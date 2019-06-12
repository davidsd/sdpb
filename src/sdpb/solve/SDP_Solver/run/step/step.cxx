#include "../../../SDP_Solver.hxx"
#include "../../../../../Timers.hxx"

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const Block_Diagonal_Matrix &A_X_inv, const Block_Diagonal_Matrix &A_Y,
  const El::Grid &block_grid, Block_Diagonal_Matrix &eigenvectors,
  Block_Diagonal_Matrix &L, Block_Matrix &L_inv_B,
  El::DistMatrix<El::BigFloat> &Q_cholesky, Timers &timers);

void compute_search_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &eigenvectors, const Block_Diagonal_Matrix &L,
  const Block_Matrix &L_inv_B, const Block_Diagonal_Matrix &X_cholesky,
  const El::BigFloat beta, const El::BigFloat &mu,
  const Block_Vector &primal_residue_p, const bool &is_corrector_phase,
  const El::DistMatrix<El::BigFloat> &Q_cholesky, Block_Vector &dx,
  Block_Diagonal_Matrix &dX, Block_Vector &dy, Block_Diagonal_Matrix &dY);

El::BigFloat
predictor_centering_parameter(const SDP_Solver_Parameters &parameters,
                              const bool is_primal_dual_feasible);

El::BigFloat corrector_centering_parameter(
  const SDP_Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_num_rows);

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma,
            const std::string &timer_name, Timers &timers);

void SDP_Solver::step(
  const SDP_Solver_Parameters &parameters, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const El::Grid &grid,
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y_cholesky,
  const Block_Diagonal_Matrix &A_X_inv, const Block_Diagonal_Matrix &A_Y,
  const Block_Vector &primal_residue_p, El::BigFloat &mu,
  El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
  El::BigFloat &dual_step_length, bool &terminate_now, Timers &timers)
{
  auto &step_timer(timers.add_and_start("run.step"));
  El::BigFloat beta_predictor;

  // Search direction: These quantities have the same structure
  // as (x, X, y, Y). They are computed twice each iteration:
  // once in the predictor step, and once in the corrector step.
  Block_Vector dx(x), dy(y);
  Block_Diagonal_Matrix dX(X), dY(Y);
  {
    // L = Cholesky decomposition of S.
    Block_Diagonal_Matrix L(block_info.schur_block_sizes,
                            block_info.block_indices,
                            block_info.schur_block_sizes.size(), grid);
    Block_Diagonal_Matrix eigenvectors(L);
    Block_Matrix L_inv_B;

    // Q = B L^-T L^-1 B
    //   = Q_cholesky Q_cholesky^T
    El::DistMatrix<El::BigFloat> Q_cholesky(sdp.dual_objective_b.Height(),
                                            sdp.dual_objective_b.Height());

    initialize_schur_complement_solver(block_info, sdp, A_X_inv, A_Y, grid,
                                       eigenvectors, L, L_inv_B, Q_cholesky,
                                       timers);

    // Compute the complementarity mu = Tr(X Y)/X.dim
    auto &frobenius_timer(
      timers.add_and_start("run.step.frobenius_product_symmetric"));
    mu = frobenius_product_symmetric(X, Y) / total_psd_rows;
    frobenius_timer.stop();
    if(mu > parameters.max_complementarity)
      {
        terminate_now = true;
        return;
      }

    auto &predictor_timer(
      timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));

    // Compute the predictor solution for (dx, dX, dy, dY)
    beta_predictor
      = predictor_centering_parameter(parameters, is_primal_and_dual_feasible);
    compute_search_direction(block_info, sdp, *this, eigenvectors, L, L_inv_B,
                             X_cholesky, beta_predictor, mu, primal_residue_p,
                             false, Q_cholesky, dx, dX, dy, dY);
    predictor_timer.stop();

    // Compute the corrector solution for (dx, dX, dy, dY)
    auto &corrector_timer(
      timers.add_and_start("run.step.computeSearchDirection(betaCorrector)"));
    beta_corrector = corrector_centering_parameter(
      parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
      total_psd_rows);

    compute_search_direction(block_info, sdp, *this, eigenvectors, L, L_inv_B,
                             X_cholesky, beta_corrector, mu, primal_residue_p,
                             true, Q_cholesky, dx, dX, dy, dY);
    corrector_timer.stop();
  }
  // Compute step-lengths that preserve positive definiteness of X, Y
  primal_step_length
    = step_length(X_cholesky, dX, parameters.step_length_reduction,
                  "run.step.stepLength(XCholesky)", timers);

  dual_step_length
    = step_length(Y_cholesky, dY, parameters.step_length_reduction,
                  "run.step.stepLength(YCholesky)", timers);

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
  step_timer.stop();
}
