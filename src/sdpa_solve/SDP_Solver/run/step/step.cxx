#include "update_cond_numbers.hxx"
#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/SDP_Solver/run/step/centering_parameters.hxx"
#include "sdp_solve/SDP_Solver/run/step/compute_R_error.hxx"
#include "sdp_solve/SDP_Solver/run/step/step_length.hxx"
#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpa_solve/memory_estimates.hxx"

namespace Sdpb::Sdpa
{
  void compute_S(const Environment &env, const SDP &sdp,
                 const Block_Info &block_info,
                 const Block_Diagonal_Matrix &X_cholesky,
                 const Block_Diagonal_Matrix &Y_cholesky, const El::Grid &grid,
                 BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
                 El::DistMatrix<El::BigFloat> &S, Timers &timers,
                 El::Matrix<int32_t> &block_timings_ms, Verbosity verbosity);

  void
  compute_search_direction(const SDP &sdp, const SDP_Solver &solver,
                           const Block_Diagonal_Matrix &minus_XY,
                           const Block_Diagonal_Matrix &X_cholesky,
                           const El::BigFloat &beta, const El::BigFloat &mu,
                           const bool &is_corrector_phase,
                           const El::DistMatrix<El::BigFloat> &S,
                           Primal_Dist_Vector &dx, Block_Diagonal_Matrix &dX,
                           Block_Diagonal_Matrix &dY);

  void SDP_Solver::step(
    const Environment &env, const Solver_Parameters &parameters,
    const Verbosity &verbosity, const std::size_t &total_psd_rows,
    const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
    const SDP &sdp, const El::Grid &grid,
    const Block_Diagonal_Matrix &X_cholesky,
    const Block_Diagonal_Matrix &Y_cholesky,
    BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context, El::BigFloat &mu,
    El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
    El::BigFloat &dual_step_length, bool &terminate_now, Timers &timers,
    El::Matrix<int32_t> &block_timings_ms, El::BigFloat &S_cond_number,
    El::BigFloat &max_block_cond_number,
    std::string &max_block_cond_number_name)
  {
    Scoped_Timer step_timer(timers, "step");
    block_timings_ms.Resize(block_info.block_dimensions.size(), 1);
    El::Zero(block_timings_ms);

    El::BigFloat beta_predictor;

    // Search direction: These quantities have the same structure
    // as (x, X, y, Y). They are computed twice each iteration:
    // once in the predictor step, and once in the corrector step.
    Primal_Dist_Vector dx(x);
    Block_Diagonal_Matrix dX(X), dY(Y);
    if(verbosity >= Verbosity::trace)
      {
        print_allocation_message_per_node(env, "dx", get_allocated_bytes(dx));
        print_allocation_message_per_node(env, "dX", get_allocated_bytes(dX));
        print_allocation_message_per_node(env, "dY", get_allocated_bytes(dY));
      }
    {
      // S_ij = Tr(G_i^T G_j), where
      // G_i = L_x_inv F_i L_Y
      // L_X, L_Y are Cholesky decomposition matrices for X and Y, respectively
      // S is an (m x m) matrix.
      El::DistMatrix<El::BigFloat> S(sdp.primal_dimension(),
                                     sdp.primal_dimension());
      if(verbosity >= Verbosity::trace)
        {
          print_allocation_message_per_node(env, "S", get_allocated_bytes(S));
        }

      compute_S(env, sdp, block_info, X_cholesky, Y_cholesky, grid,
                bigint_syrk_context, S, timers, block_timings_ms, verbosity);
      try
        {
          Scoped_Timer Cholesky_timer(timers, "Cholesky_S");
          Cholesky(El::UpperOrLowerNS::UPPER, S);
        }
      catch(std::exception &e)
        {
          RUNTIME_ERROR("Error when computing Cholesky(S): ", e.what());
        }

      // Calculate matrix product -XY
      // It will be reused for mu, R-err, compute_search_direction().
      Scoped_Timer XY_timer(timers, "XY");
      Block_Diagonal_Matrix minus_XY(X);
      if(verbosity >= Verbosity::trace)
        {
          print_allocation_message_per_node(env, "XY",
                                            get_allocated_bytes(minus_XY));
        }
      scale_multiply_add(El::BigFloat(-1), X, Y, El::BigFloat(0), minus_XY);
      XY_timer.stop();

      // Compute the complementarity mu = Tr(X Y)/X.dim
      {
        Scoped_Timer mu_timer(timers, "mu");
        mu = -minus_XY.trace() / total_psd_rows;
      }
      if(mu > parameters.max_complementarity)
        {
          terminate_now = true;
          // Block timings are not updated below, so here we already have correct values.
          // (otherwise, we might want to clear them)
          Scoped_Timer block_timings_timer(timers, "block_timings_AllReduce");
          El::AllReduce(block_timings_ms, El::mpi::COMM_WORLD);
          return;
        }

      // R = mu * I - XY
      // R_error = maxAbs(R)
      R_error = compute_R_error(mu, minus_XY, timers);

      {
        Scoped_Timer predictor_timer(timers,
                                     "computeSearchDirection(betaPredictor)");

        // Compute the predictor solution for (dx, dX, dy, dY)
        beta_predictor = predictor_centering_parameter(
          parameters, is_primal_and_dual_feasible);
        compute_search_direction(sdp, *this, minus_XY, X_cholesky,
                                 beta_predictor, mu, false, S, dx, dX, dY);
      }

      // Compute the corrector solution for (dx, dX, dy, dY)
      {
        Scoped_Timer corrector_timer(timers,
                                     "computeSearchDirection(betaCorrector)");
        beta_corrector = corrector_centering_parameter(
          parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
          total_psd_rows);

        compute_search_direction(sdp, *this, minus_XY, X_cholesky,
                                 beta_corrector, mu, true, S, dx, dX, dY);
      }

      // Calculate condition numbers for Cholesky matrices
      update_cond_numbers(S, block_info, X_cholesky, Y_cholesky, timers,
                          S_cond_number, max_block_cond_number,
                          max_block_cond_number_name);
    }
    // Compute step-lengths that preserve positive definiteness of X, Y
    primal_step_length
      = step_length(X_cholesky, dX, parameters.step_length_reduction,
                    "stepLength(XCholesky)", timers);

    dual_step_length
      = step_length(Y_cholesky, dY, parameters.step_length_reduction,
                    "stepLength(YCholesky)", timers);

    // If our problem is both dual-feasible and primal-feasible,
    // ensure we're following the true Newton direction.
    if(is_primal_and_dual_feasible)
      {
        primal_step_length = El::Min(primal_step_length, dual_step_length);
        dual_step_length = primal_step_length;
      }

    // Update the primal point (x, X) += primalStepLength*(dx, dX)
    El::Axpy(primal_step_length, dx, x);
    dX *= primal_step_length;
    X += dX;

    // Update the dual point Y += dualStepLength * dY
    dY *= dual_step_length;
    Y += dY;

    // Block timings
    Scoped_Timer block_timings_timer(timers, "block_timings_AllReduce");
    El::AllReduce(block_timings_ms, El::mpi::COMM_WORLD);
  }
}
