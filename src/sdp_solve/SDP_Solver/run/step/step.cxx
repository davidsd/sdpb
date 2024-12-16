#include "compute_R_error.hxx"
#include "update_cond_numbers.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpb_util/memory_estimates.hxx"

void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C);
// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

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
            const std::string &timer_name, El::BigFloat &max_step,
            Timers &timers);

void compute_errors(const std::size_t &total_psd_rows,
                    const Block_Vector &x_const, const Block_Vector &dx_const,
                    const Block_Vector &y_const, const Block_Vector &dy_const,
                    const Block_Diagonal_Matrix &X_const,
                    const Block_Diagonal_Matrix &dX_const,
                    const Block_Diagonal_Matrix &Y_const,
                    const Block_Diagonal_Matrix &dY_const,
                    const El::BigFloat &primal_step_length,
                    const El::BigFloat &dual_step_length,
                    El::BigFloat &R_error, El::BigFloat &mu, Timers &timers)
{
  Scoped_Timer timer(timers, "compute_errors");
  Block_Vector x(x_const), y(y_const), dx(dx_const), dy(dy_const);
  Block_Diagonal_Matrix X(X_const), Y(Y_const), dX(dX_const), dY(dY_const);

  // Update x, y, dX, dY ///////////
  for(size_t block = 0; block < x.blocks.size(); ++block)
    {
      El::Axpy(primal_step_length, dx.blocks[block], x.blocks[block]);
    }
  dX *= primal_step_length;
  X += dX;
  for(size_t block = 0; block < dy.blocks.size(); ++block)
    {
      El::Axpy(dual_step_length, dy.blocks[block], y.blocks[block]);
    }
  dY *= dual_step_length;
  Y += dY;
  //////////////////////////////////

  mu = frobenius_product_symmetric(X, Y) / total_psd_rows;

  Block_Diagonal_Matrix R(X);
  scale_multiply_add(El::BigFloat(-1), X, Y, El::BigFloat(0), R);
  R.add_diagonal(mu);

  R_error = R.max_abs();
}

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
  El::BigFloat &dual_step_length, size_t &num_corrector_iterations,
  bool &terminate_now, Timers &timers, El::Matrix<int32_t> &block_timings_ms,
  El::BigFloat &Q_cond_number, El::BigFloat &max_block_cond_number,
  std::string &max_block_cond_number_name)
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

#define VERBOSE_ALLOCATION_MESSAGE(var)                                       \
  if(verbosity >= Verbosity::trace)                                           \
    print_allocation_message_per_node(env, #var, get_allocated_bytes(var));

  VERBOSE_ALLOCATION_MESSAGE(dx);
  VERBOSE_ALLOCATION_MESSAGE(dy);
  VERBOSE_ALLOCATION_MESSAGE(dX);
  VERBOSE_ALLOCATION_MESSAGE(dY);

  {
    // SchurComplementCholesky = L', the Cholesky decomposition of the
    // Schur complement matrix S.
    Block_Diagonal_Matrix schur_complement_cholesky(
      block_info.schur_block_sizes(), block_info.block_indices,
      block_info.num_points.size(), grid);

    VERBOSE_ALLOCATION_MESSAGE(schur_complement_cholesky);

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
    VERBOSE_ALLOCATION_MESSAGE(Q);

    // Compute SchurComplement and prepare to solve the Schur
    // complement equation for dx, dy
    initialize_schur_complement_solver(env, block_info, sdp, A_X_inv, A_Y,
                                       grid, schur_complement_cholesky,
                                       schur_off_diagonal, bigint_syrk_context,
                                       Q, timers, block_timings_ms, verbosity);

    // Calculate condition numbers for Cholesky matrices
    update_cond_numbers(Q, block_info, schur_complement_cholesky, X_cholesky,
                        Y_cholesky, timers, Q_cond_number,
                        max_block_cond_number, max_block_cond_number_name);

    // Calculate matrix product -XY
    // It will be reused for mu, R-err, compute_search_direction().
    Scoped_Timer XY_timer(timers, "XY");
    Block_Diagonal_Matrix minus_XY(X);
    VERBOSE_ALLOCATION_MESSAGE(minus_XY);

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
    //
    // TODO: now we always update R_error during corrector phase by calling compute_errors:
    // R = mu' * I - X'Y'
    //   where:
    //   X' = X + dX, Y' = Y + dY
    //   mu' = Tr(X'Y') / X'.dim
    // TODO: which definition should we use?
    R_error = compute_R_error(mu, minus_XY, timers);

    // If set to 'true', then we perform a centering step,
    // i.e. with fixed mu (beta = 1)
    // TODO: now we never perform a centering step.
    // In principle, we could do it if we are far from the local central path
    // (i.e. R-err is large)
    const bool do_centering_step = false;
    // The condition below leads to solver stalling at mu ~ 1e37
    //   for SingletScalar_cT_test_nmax6 test case:
    // const bool do_centering_step = R_error > mu;

    {
      Scoped_Timer predictor_timer(timers, "predictor");

      // Compute the predictor solution for (dx, dX, dy, dY)
      beta_predictor = predictor_centering_parameter(
        parameters, is_primal_and_dual_feasible);

      if(do_centering_step)
        beta_predictor = 1;

      compute_search_direction(block_info, sdp, *this, minus_XY,
                               schur_complement_cholesky, schur_off_diagonal,
                               X_cholesky, beta_predictor, mu,
                               primal_residue_p, false, Q, dx, dX, dy, dY);
    }

    // Compute the corrector solution for (dx, dX, dy, dY)
    {
      Scoped_Timer corrector_timer(timers, "corrector");

      Block_Vector dx_prev(dx), dy_prev(dy);
      Block_Diagonal_Matrix dX_prev(dX), dY_prev(dY);

      VERBOSE_ALLOCATION_MESSAGE(dx_prev);
      VERBOSE_ALLOCATION_MESSAGE(dy_prev);
      VERBOSE_ALLOCATION_MESSAGE(dX_prev);
      VERBOSE_ALLOCATION_MESSAGE(dY_prev);

      // Will be initialized at the end of the first corrector iteration
      El::BigFloat primal_step_length_prev = 0;
      El::BigFloat dual_step_length_prev = 0;
      El::BigFloat beta_corrector_prev = 0;

      El::BigFloat reduce_factor
        = 1
          - El::Min(primal_step_length, dual_step_length)
              * (1 - beta_corrector);
      El::BigFloat reduce_factor_prev = 1;

      beta_corrector = corrector_centering_parameter(
        parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
        total_psd_rows);

      const El::BigFloat corrector_iter_mu_reduction
        = parameters.corrector_mu_reduction;
      ASSERT(corrector_iter_mu_reduction > 0,
             DEBUG_STRING(corrector_iter_mu_reduction));
      ASSERT(corrector_iter_mu_reduction < 1,
             DEBUG_STRING(corrector_iter_mu_reduction));

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

      bool undo_last_corrector_iteration = false;
      num_corrector_iterations = 0;
      while(num_corrector_iterations < max_corrector_iterations)
        {
          Scoped_Timer loop_timer(timers,
                                  std::to_string(num_corrector_iterations));

          if(num_corrector_iterations > 0)
            {
              dx_prev = dx;
              dy_prev = dy;
              dX_prev = dX;
              dY_prev = dY;
              primal_step_length_prev = primal_step_length;
              dual_step_length_prev = dual_step_length;
              beta_corrector_prev = beta_corrector;
              reduce_factor_prev = reduce_factor;

              beta_corrector = beta_corrector * corrector_iter_mu_reduction;
            }

          Scoped_Timer compute_search_direction_timer(
            timers, "compute_search_direction");
          compute_search_direction(
            block_info, sdp, *this, minus_XY, schur_complement_cholesky,
            schur_off_diagonal, X_cholesky, beta_corrector, mu,
            primal_residue_p, true, Q, dx, dX, dy, dY);
          compute_search_direction_timer.stop();

          // Compute step-lengths that preserve positive definiteness of X, Y
          El::BigFloat max_primal_step_length;
          primal_step_length = step_length(
            X_cholesky, dX, parameters.step_length_reduction,
            "stepLength(XCholesky)", max_primal_step_length, timers);
          El::BigFloat max_dual_step_length;
          dual_step_length = step_length(
            Y_cholesky, dY, parameters.step_length_reduction,
            "stepLength(YCholesky)", max_dual_step_length, timers);

          // Update R-err,
          // print corrector iteration status
          {
            El::BigFloat error_P, error_p, error_d, coit_mu;

            compute_errors(total_psd_rows, x, dx, y, dy, X, dX, Y, dY,
                           primal_step_length, dual_step_length, this->R_error,
                           coit_mu, timers);

            const auto min_step_length
              = El::Min(primal_step_length, dual_step_length);
            reduce_factor = 1 - min_step_length * (1 - beta_corrector);
            if(El::mpi::Rank() == 0 && verbosity >= Verbosity::debug)
              {
                El::Output("  step=(", primal_step_length, ",",
                           dual_step_length, ") maxstep=(",
                           max_primal_step_length, ",", max_dual_step_length,
                           ") R=", R_error, " mu=", coit_mu,
                           " reduce=", reduce_factor);
              }
          }
          ++num_corrector_iterations;

          // continue corrector steps
          // only if (primal_step_length + dual_step_length)
          // is not too small compared to the historical maximum.
          // TODO: check other exit conditions, e.g (gap < dualityGapThreshold)

          if(El::Max(primal_step_length, dual_step_length)
             < parameters.corrector_step_length_threshold)
            {
              if(num_corrector_iterations > 1)
                undo_last_corrector_iteration = true;
              // TODO in this case, we can exit even before computing and printing errors, shall we?
              break;
            }

          if(num_corrector_iterations > 1
             && reduce_factor >= reduce_factor_prev)
            {
              undo_last_corrector_iteration = true;
              break;
            }
        }

      if(El::mpi::Rank() == 0 && verbosity >= Verbosity::debug)
        {
          El::Output(
            "  num_corrector_iterations=", num_corrector_iterations,
            undo_last_corrector_iteration ? ", the last one discarded." : "");
        }

      if(undo_last_corrector_iteration)
        {
          ASSERT(num_corrector_iterations > 0,
                 "Should perform at least one corrector iteration!");
          dx = dx_prev;
          dy = dy_prev;
          dX = dX_prev;
          dY = dY_prev;
          primal_step_length = primal_step_length_prev;
          dual_step_length = dual_step_length_prev;
          beta_corrector = beta_corrector_prev;
          reduce_factor = reduce_factor_prev;
        }
    }
  }

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

#undef VERBOSE_ALLOCATION_MESSAGE
}
