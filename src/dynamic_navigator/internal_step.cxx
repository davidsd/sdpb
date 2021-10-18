#include "../../../SDP_Solver.hxx"


El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

El::BigFloat predictor_centering_parameter(const Solver_Parameters &parameters,
                                           const bool is_primal_dual_feasible);

void internal_predictor_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);


//Compute the predictor step of the center sdp 
//Return H_xx^-1 grad_x(Lag)
void internal_step(
  const Solver_Parameters &parameters, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal, const El::DistMatrix<El::BigFloat> &Q,
  const Block_Vector &primal_residue_p, El::BigFloat &mu,Timers &timers,
  Block_Vector &dx, Block_Vector &dy)
{
  auto &step_timer(timers.add_and_start("run.step"));
  El::BigFloat beta_predictor;

  // Search direction: These quantities have the same structure
  // as (x, X, y, Y). They are computed twice each iteration:
  // once in the predictor step, and once in the corrector step.
  Block_Vector dx(x), dy(y);
  Block_Diagonal_Matrix dX(X), dY(Y);
  
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


  auto &cholesky_decomposition_timer(
      timers.add_and_start("run.choleskyDecomposition"));
  Block_Diagonal_Matrix X_cholesky(X)
  cholesky_decomposition(X, X_cholesky);


  // Solve the predictor solution for (dx,dy)
  auto &predictor_timer(
    timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));
  beta_predictor
    = predictor_centering_parameter(parameters, is_primal_and_dual_feasible);

  internal_predictor_direction(block_info, sdp, *this, schur_complement_cholesky,
                           schur_off_diagonal, X_cholesky, beta_predictor,
                           mu, primal_residue_p, Q, dx, dy);
}




