#include "Corrector_Iteration.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/memory_estimates.hxx"

#define VERBOSE_ALLOCATION_MESSAGE(var)                                       \
  if(verbosity >= Verbosity::trace)                                           \
    print_allocation_message_per_node(env, #var, get_allocated_bytes(var));

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C);

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

El::BigFloat corrector_centering_parameter(
  const Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_num_rows);

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma,
            const El::BigFloat &boost_step_min,
            const El::BigFloat &boost_step_max, const std::string &timer_name,
            El::BigFloat &max_step, Timers &timers);

// R = mu * I - X' Y',
// where
//   mu = Tr(X'Y') / dim(X')
//   X' = X + primal_step_length * dX
//   Y' = Y + primal_step_length * dY
Block_Diagonal_Matrix
compute_R_matrix(const std::size_t &total_psd_rows, const SDP_Solver &solver,
                 const Block_Diagonal_Matrix &dX_const,
                 const Block_Diagonal_Matrix &dY_const,
                 const El::BigFloat &primal_step_length,
                 const El::BigFloat &dual_step_length, El::BigFloat &mu)
{
  // X -> X + P_step * dX
  Block_Diagonal_Matrix X = solver.X;
  for(size_t i = 0; i < X.blocks.size(); ++i)
    {
      El::Axpy(primal_step_length, dX_const.blocks.at(i), X.blocks.at(i));
    }

  // Y -> Y + D_step * dY
  Block_Diagonal_Matrix Y = solver.Y;
  for(size_t i = 0; i < Y.blocks.size(); ++i)
    {
      El::Axpy(dual_step_length, dY_const.blocks.at(i), Y.blocks.at(i));
    }

  // mu = Tr(XY) / dim(X)
  mu = frobenius_product_symmetric(X, Y) / total_psd_rows;

  // R = mu * I - XY,
  Block_Diagonal_Matrix R(X);
  scale_multiply_add(El::BigFloat(-1), X, Y, El::BigFloat(0), R);
  R.add_diagonal(mu);
  return R;
}

Corrector_Iteration single_corrector_iteration(
  const SDP_Solver &solver, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q,
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y_cholesky,
  const Block_Diagonal_Matrix &minus_XY, const Block_Vector &primal_residue_p,
  const El::BigFloat &beta_corrector, const Solver_Parameters &parameters,
  // mu = Tr(X Y)/X.dim
  const El::BigFloat &initial_mu, Block_Vector &dx, Block_Diagonal_Matrix &dX,
  Block_Vector &dy, Block_Diagonal_Matrix &dY, Timers &timers)
{
  Corrector_Iteration iter;

  {
    Scoped_Timer compute_search_direction_timer(timers,
                                                "compute_search_direction");
    constexpr bool is_corrector_phase = true;
    compute_search_direction(
      block_info, sdp, solver, minus_XY, schur_complement_cholesky,
      schur_off_diagonal, X_cholesky, beta_corrector, initial_mu,
      primal_residue_p, is_corrector_phase, Q, dx, dX, dy, dY);
  }

  {
    iter.primal_step_length = step_length(
      X_cholesky, dX, parameters.step_length_reduction,
      parameters.corrector_step_boost_min, parameters.corrector_step_boost_max,
      "stepLength(XCholesky)", iter.max_primal_step_length, timers);
    iter.dual_step_length = step_length(
      Y_cholesky, dY, parameters.step_length_reduction,
      parameters.corrector_step_boost_min, parameters.corrector_step_boost_max,
      "stepLength(YCholesky)", iter.max_dual_step_length, timers);
    // If our problem is both dual-feasible and primal-feasible,
    // ensure we're following the true Newton direction.
    if(is_primal_and_dual_feasible)
      {
        iter.primal_step_length = iter.dual_step_length
          = El::Min(iter.primal_step_length, iter.dual_step_length);
      }
  }

  {
    Scoped_Timer timer(timers, "compute_errors");
    const auto R = compute_R_matrix(total_psd_rows, solver, dX, dY,
                                    iter.primal_step_length,
                                    iter.dual_step_length, iter.mu);

    iter.R_error = R.max_abs();
    iter.R_mean_abs = R.mean_abs();

    {
      Scoped_Timer dXdY_timer(timers, "dXdY");
      Block_Diagonal_Matrix dXdY(dX);
      scale_multiply_add(1, dX, dY, 0, dXdY);
      El::BigFloat total = 0;
      int size = 0;
      for(const auto &block : dXdY.blocks)
        {
          const auto &m = block.LockedMatrix();
          size += m.Height() * m.Width();
          for(int i = 0; i < m.Height(); ++i)
            for(int j = 0; j < m.Width(); ++j)
              {
                // TODO should we also multiply by P-step and D-step?
                total += m(i, j);
              };
        }
      total = El::mpi::AllReduce(total, El::mpi::SUM, El::mpi::COMM_WORLD);
      size = El::mpi::AllReduce(size, El::mpi::SUM, El::mpi::COMM_WORLD);
      iter.dXdY_mean = total / size;
    }
  }

  return iter;
}

// -d(log10 mu)/dt, with t measured in ms
// Optionally multiplied by normalization factor
Boost_Float
get_log_mu_speed(const El::BigFloat &mu_0, const El::BigFloat &mu_1,
                 int64_t time_ms, const Boost_Float &normalization)
{
  if(time_ms == 0)
    time_ms = 1;

  ASSERT(time_ms > 0, DEBUG_STRING(time_ms));
  return log10(to_Boost_Float(mu_0) / to_Boost_Float(mu_1)) / time_ms
         * normalization;
}

void corrector_step(
  const SDP_Solver &solver, const Scoped_Timer &step_timer,
  const Environment &env, const Solver_Parameters &parameters,
  const Verbosity &verbosity, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q,
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y_cholesky,
  const Block_Diagonal_Matrix &minus_XY, const Block_Vector &primal_residue_p,
  const bool do_centering_step, const El::BigFloat &mu,
  El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
  El::BigFloat &dual_step_length, Block_Vector &dx, Block_Vector &dy,
  Block_Diagonal_Matrix &dX, Block_Diagonal_Matrix &dY,
  std::vector<Corrector_Iteration> &corrector_iterations, Timers &timers)
{
  Scoped_Timer corrector_timer(timers, "corrector");

  Block_Vector dx_prev(dx), dy_prev(dy);
  Block_Diagonal_Matrix dX_prev(dX), dY_prev(dY);

  VERBOSE_ALLOCATION_MESSAGE(dx_prev);
  VERBOSE_ALLOCATION_MESSAGE(dy_prev);
  VERBOSE_ALLOCATION_MESSAGE(dX_prev);
  VERBOSE_ALLOCATION_MESSAGE(dY_prev);

  beta_corrector = corrector_centering_parameter(
    parameters, solver.X, dX, solver.Y, dY, mu, is_primal_and_dual_feasible,
    total_psd_rows);

  El::BigFloat corrector_iter_mu_reduction = parameters.corrector_mu_reduction;
  if(corrector_iter_mu_reduction <= 0)
    corrector_iter_mu_reduction = beta_corrector;

  size_t max_corrector_iterations
    = is_primal_and_dual_feasible
        ? parameters.feasible_max_corrector_iterations
        : parameters.infeasible_max_corrector_iterations;
  if(max_corrector_iterations == 0)
    max_corrector_iterations = std::numeric_limits<int64_t>::max();

  if(do_centering_step)
    {
      beta_corrector = 1;
      max_corrector_iterations = 1;
    }

  // How fast do we decrease mu by doing one predictor + one corrector step
  Boost_Float log_mu_speed_normalization;

  bool undo_last_corrector_iteration = false;
  corrector_iterations.clear();
  while(corrector_iterations.size() < max_corrector_iterations)
    {
      Scoped_Timer loop_timer(timers,
                              std::to_string(corrector_iterations.size()));

      if(!corrector_iterations.empty())
        {
          dx_prev = dx;
          dy_prev = dy;
          dX_prev = dX;
          dY_prev = dY;

          // Aim for lower mu each time
          beta_corrector = beta_corrector * corrector_iter_mu_reduction;
        }

      // Compute dx,dy,dX,dY etc.
      auto &iteration
        = corrector_iterations.emplace_back(single_corrector_iteration(
          solver, total_psd_rows, is_primal_and_dual_feasible, block_info, sdp,
          schur_complement_cholesky, schur_off_diagonal, Q, X_cholesky,
          Y_cholesky, minus_XY, primal_residue_p, beta_corrector, parameters,
          mu, dx, dX, dy, dY, timers));

      // Previous iteration
      auto prev_it = corrector_iterations.rbegin();
      ++prev_it;

      const int64_t total_step_time_ms = step_timer.elapsed_milliseconds();
      if(corrector_iterations.size() == 1)
        {
          // Normalization:
          // We set speed = 1
          // If we decreased mu by a factor of 10
          // during one full solver step (1 predictor + 1 corrector)
          log_mu_speed_normalization
            = 1 / get_log_mu_speed(mu, mu / 10, total_step_time_ms, 1);
        }
      iteration.log_mu_speed_full = get_log_mu_speed(
        mu, iteration.mu, total_step_time_ms, log_mu_speed_normalization);
      {
        auto prev_mu = corrector_iterations.size() == 1 ? mu : prev_it->mu;
        iteration.log_mu_speed_corrector = get_log_mu_speed(
          prev_mu, iteration.mu, loop_timer.elapsed_milliseconds(),
          log_mu_speed_normalization);
      }

      // print corrector iteration status
      if(El::mpi::Rank() == 0 && verbosity >= Verbosity::debug)
        {
          El::Output("  step=(", iteration.primal_step_length, ",",
                     iteration.dual_step_length, ") maxstep=(",
                     iteration.max_primal_step_length, ",",
                     iteration.max_dual_step_length, ") R=", iteration.R_error,
                     " R_mean=", iteration.R_mean_abs,
                     " dXdY_mean=", iteration.dXdY_mean, " mu=", iteration.mu,
                     " -dlog(mu)/dt_corr=", iteration.log_mu_speed_corrector,
                     " -dlog(mu)/dt_full=", iteration.log_mu_speed_full);
        }

      // Check whether we should continue corrector iterations.
      // TODO: check also global exit conditions, e.g (gap < dualityGapThreshold)

      // Stop if mu does not decrease
      if(corrector_iterations.size() > 1 && iteration.mu >= prev_it->mu)
        {
          undo_last_corrector_iteration = true;
          break;
        }

      // Continue corrector iterations
      // only if max(primal_step_length, dual_step_length) is not too small.
      if(El::Max(iteration.primal_step_length, iteration.dual_step_length)
         < parameters.corrector_step_length_threshold)
        {
          if(corrector_iterations.size() > 1)
            undo_last_corrector_iteration = true;
          break;
        }

      // If current corrector decreases mu too slowly,
      // then we should make full solver step again.
      if(parameters.corrector_check_mu_speed && corrector_iterations.size() > 1
         && iteration.log_mu_speed_corrector
              < corrector_iterations.front().log_mu_speed_full)
        {
          break;
        }
    }

  if(El::mpi::Rank() == 0 && verbosity >= Verbosity::debug)
    {
      El::Output("  num_corrector_iterations=", corrector_iterations.size(),
                 undo_last_corrector_iteration ? ", the last one discarded."
                                               : "");
    }

  ASSERT(!corrector_iterations.empty(),
         "Expected at least one corrector iteration!");
  auto last_successful_iteration = corrector_iterations.back();
  if(undo_last_corrector_iteration)
    {
      corrector_iterations.back().is_canceled = true;
      ASSERT(corrector_iterations.size() >= 2,
             DEBUG_STRING(corrector_iterations.size()),
             "The first corrector iteration cannot be undone!");
      dx = dx_prev;
      dy = dy_prev;
      dX = dX_prev;
      dY = dY_prev;
      last_successful_iteration
        = corrector_iterations.at(corrector_iterations.size() - 2);
    }

  // Initialize output variables

  const bool extra_corrector_success = [&] {
    if(corrector_iterations.size() == 1)
      return false;
    if(corrector_iterations.size() == 2 && undo_last_corrector_iteration)
      return false;
    if(parameters.corrector_R_threshold > 0)
      {
        if(last_successful_iteration.R_mean_abs
           > corrector_iterations.front().R_mean_abs
               * parameters.corrector_R_threshold)
          {
            return false;
          }
      }
    return true;
  }();

  if(extra_corrector_success)
    {
      // Use full step size
      primal_step_length = last_successful_iteration.primal_step_length;
      dual_step_length = last_successful_iteration.dual_step_length;
    }
  else
    {
      // Use reduced step size, as in the ordinary SDPB algorithm.
      primal_step_length
        = El::Min(last_successful_iteration.max_primal_step_length
                    * parameters.step_length_reduction,
                  El::BigFloat(1));
      dual_step_length = El::Min(last_successful_iteration.max_dual_step_length
                                   * parameters.step_length_reduction,
                                 El::BigFloat(1));

      if(El::mpi::Rank() == 0 && verbosity >= Verbosity::debug)
        {
          El::Output("  reset step=(", primal_step_length, ",",
                     dual_step_length, ")");
        }
    }
}

#undef VERBOSE_ALLOCATION_MESSAGE