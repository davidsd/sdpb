#include "compute_R_error.hxx"
#include "update_cond_numbers.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpb_util/memory_estimates.hxx"

void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C);

void initialize_schur_complement_solver(
  const Environment &env, const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const El::Grid &block_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal,
  BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
  El::DistMatrix<El::BigFloat> &Q, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms, Verbosity verbosity);

void compute_search_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &minus_XY,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat &beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const bool &is_corrector_phase, const El::DistMatrix<El::BigFloat> &Q,
  Block_Vector &dx, Block_Diagonal_Matrix &dX, Block_Vector &dy,
  Block_Diagonal_Matrix &dY);

El::BigFloat predictor_centering_parameter(const Solver_Parameters &parameters,
                                           const bool is_primal_dual_feasible);

El::BigFloat corrector_centering_parameter(
  const Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_num_rows);

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma,
            const std::string &timer_name, Timers &timers);

void SDP_Solver::step(
  const Environment &env, const Solver_Parameters &parameters,
  const Verbosity &verbosity, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const El::Grid &grid,
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y_cholesky,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const Block_Vector &primal_residue_p,
  BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context, El::BigFloat &mu,
  El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
  El::BigFloat &dual_step_length, bool &terminate_now, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms, El::BigFloat &Q_cond_number,
  El::BigFloat &max_block_cond_number, std::string &max_block_cond_number_name)
{
  Scoped_Timer step_timer(timers, "step");
  block_timings_ms.Resize(block_info.dimensions.size(), 1);
  El::Zero(block_timings_ms);

  El::BigFloat beta_predictor;

  // Search direction: These quantities have the same structure
  // as (x, X, y, Y). They are computed twice each iteration:
  // once in the predictor step, and once in the corrector step.
  Block_Vector dx(x), dy(y);
  Block_Diagonal_Matrix dX(X), dY(Y);
  if(verbosity >= Verbosity::trace)
    {
      print_allocation_message_per_node(env, "dx", get_allocated_bytes(dx));
      print_allocation_message_per_node(env, "dy", get_allocated_bytes(dy));
      print_allocation_message_per_node(env, "dX", get_allocated_bytes(dX));
      print_allocation_message_per_node(env, "dY", get_allocated_bytes(dY));
    }
  {
    // SchurComplementCholesky = L', the Cholesky decomposition of the
    // Schur complement matrix S.
    Block_Diagonal_Matrix schur_complement_cholesky(
      block_info.schur_block_sizes(), block_info.block_indices,
      block_info.num_points.size(), grid);
    if(verbosity >= Verbosity::trace)
      {
        print_allocation_message_per_node(
          env, "schur_complement_cholesky",
          get_allocated_bytes(schur_complement_cholesky));
      }

    // SchurOffDiagonal = L'^{-1} FreeVarMatrix, needed in solving the
    // Schur complement equation.
    Block_Matrix schur_off_diagonal;

    // Q = B' L'^{-T} L'^{-1} B' - {{0, 0}, {0, 1}}, where B' =
    // (FreeVarMatrix U).  Q is needed in the factorization of the Schur
    // complement equation.  Q has dimension N'xN', where
    //
    //   N' = cols(B) + cols(U) = N + cols(U)
    //
    // where N is the dimension of the dual objective function.  Note
    // that N' could change with each iteration.
    El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
                                   sdp.dual_objective_b.Height());
    if(verbosity >= Verbosity::trace)
      {
        print_allocation_message_per_node(env, "Q", get_allocated_bytes(Q));
      }

    // Compute SchurComplement and prepare to solve the Schur
    // complement equation for dx, dy
    initialize_schur_complement_solver(env, block_info, sdp, A_X_inv, A_Y,
                                       grid, schur_complement_cholesky,
                                       schur_off_diagonal, bigint_syrk_context,
                                       Q, timers, block_timings_ms, verbosity);

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
      compute_search_direction(block_info, sdp, *this, minus_XY,
                               schur_complement_cholesky, schur_off_diagonal,
                               X_cholesky, beta_predictor, mu,
                               primal_residue_p, false, Q, dx, dX, dy, dY);
    }

    // Compute the corrector solution for (dx, dX, dy, dY)
    {
      Scoped_Timer corrector_timer(timers,
                                   "computeSearchDirection(betaCorrector)");
      beta_corrector = corrector_centering_parameter(
        parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
        total_psd_rows);

      compute_search_direction(block_info, sdp, *this, minus_XY,
                               schur_complement_cholesky, schur_off_diagonal,
                               X_cholesky, beta_corrector, mu,
                               primal_residue_p, true, Q, dx, dX, dy, dY);
    }

    // Calculate condition numbers for Cholesky matrices
    update_cond_numbers(Q, block_info, schur_complement_cholesky, X_cholesky,
                        Y_cholesky, timers, Q_cond_number,
                        max_block_cond_number, max_block_cond_number_name);
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
  for(size_t block = 0; block < x.blocks.size(); ++block)
    {
      El::Axpy(primal_step_length, dx.blocks[block], x.blocks[block]);
    }
  dX *= primal_step_length;

  X += dX;

  // Update the dual point (y, Y) += dualStepLength*(dy, dY)
  for(size_t block = 0; block < dy.blocks.size(); ++block)
    {
      El::Axpy(dual_step_length, dy.blocks[block], y.blocks[block]);
    }
  dY *= dual_step_length;

  Y += dY;

  // Block timings
  Scoped_Timer block_timings_timer(timers, "block_timings_AllReduce");
  El::AllReduce(block_timings_ms, El::mpi::COMM_WORLD);
}
