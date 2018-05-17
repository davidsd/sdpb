#include "../constraint_matrix_weighted_sum.hxx"
#include "../../../../Timers.hxx"

// Compute the search direction (dx, dX, dy, dY) for the predictor and
// corrector phases.
//
// Inputs:
// - beta, the centering parameter
// - mu = Tr(X Y) / X.cols
// - correctorPhase: boolean indicating whether we're in the corrector
//   phase or predictor phase.
// Workspace (members of SDPSolver which are modified in-place but not
// used elsewhere):
// - Z, R
// Outputs (members of SDPSolver which are modified in-place):
// - dx, dX, dy, dY
//

void compute_schur_RHS(const SDP &sdp, const Vector &dual_residues,
                       const Block_Vector &dual_residues_elemental,
                       const Block_Diagonal_Matrix &Z, const Vector &x,
                       const Block_Vector &x_elemental, Vector &r_x,
                       Block_Vector &r_x_elemental, Vector &r_y,
                       El::DistMatrix<El::BigFloat> &r_y_elemental);

void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx,
  El::DistMatrix<El::BigFloat> &dy);

void SDP_Solver::compute_search_direction(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal_elemental,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta_elemental,
  const El::BigFloat &mu_elemental, const bool correctorPhase,
  const El::DistMatrix<El::BigFloat> &Q_elemental, Vector &dx,
  Block_Vector &dx_elemental, Block_Diagonal_Matrix &dX, Vector &dy,
  El::DistMatrix<El::BigFloat> &dy_elemental, Block_Diagonal_Matrix &dY)
{
  std::string timerName = "computeSearchDirection(";
  if(correctorPhase)
    {
      timerName += "betaCorrector)";
    }
  else
    {
      timerName += "betaPredictor)";
    }

  // R = beta mu I - X Y (predictor phase)
  // R = beta mu I - X Y - dX dY (corrector phase)
  Block_Diagonal_Matrix R(X);

  timers[timerName + ".R.XY"].resume();
  block_diagonal_matrix_scale_multiply_add(Real(-1), X, Y, Real(0), R);
  block_diagonal_matrix_scale_multiply_add(El::BigFloat(-1), X, Y,
                                           El::BigFloat(0), R);
  timers[timerName + ".R.XY"].stop();
  if(correctorPhase)
    {
      timers[timerName + ".R.dXdY"].resume();
      block_diagonal_matrix_scale_multiply_add(Real(-1), dX, dY, Real(1), R);
      block_diagonal_matrix_scale_multiply_add(El::BigFloat(-1), dX, dY,
                                               El::BigFloat(1), R);
      timers[timerName + ".R.dXdY"].stop();
    }

  timers[timerName + ".R.add"].resume();
  R.add_diagonal(beta_elemental * mu_elemental);
  timers[timerName + ".R.add"].stop();

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))

  // FIXME: It has to have the same shape as X, but does not need to
  // be a copy.
  Block_Diagonal_Matrix Z(X);

  timers[timerName + ".Z.multiply"].resume();
  block_diagonal_matrix_multiply(primal_residues, Y, Z);
  timers[timerName + ".Z.multiply"].stop();
  timers[timerName + ".Z.subtract"].resume();
  Z -= R;
  timers[timerName + ".Z.subtract"].stop();
  timers[timerName + ".Z.cholesky"].resume();
  block_matrix_solve_with_cholesky(X_cholesky, Z);
  timers[timerName + ".Z.cholesky"].stop();
  timers[timerName + ".Z.symm"].resume();
  Z.symmetrize();
  timers[timerName + ".Z.symm"].stop();

  // r_x[p] = -dual_residues[p] - Tr(A_p Z)
  // r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
  // Here, dx = r_x, dy = r_y.
  timers[timerName + ".computeSchurRHS"].resume();
  compute_schur_RHS(sdp, dual_residues, dual_residues_elemental, Z, x,
                    x_elemental, dx, dx_elemental, dy, dy_elemental);
  timers[timerName + ".computeSchurRHS"].stop();

  // Solve for dx, dy in-place
  timers[timerName + ".dxdy"].resume();
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal_elemental, Q_elemental,
                                  dx_elemental, dy_elemental);
  timers[timerName + ".dxdy"].stop();

  // dX = PrimalResidues + \sum_p A_p dx[p]
  timers[timerName + ".dX.weightedSum"].resume();
  constraint_matrix_weighted_sum(sdp, dx, dX);
  constraint_matrix_weighted_sum(sdp, dx_elemental, dX);
  timers[timerName + ".dX.weightedSum"].stop();
  timers[timerName + ".dX.primalRes"].resume();
  dX += primal_residues;
  timers[timerName + ".dX.primalRes"].stop();

  // dY = Symmetrize(X^{-1} (R - dX Y))
  timers[timerName + ".dY.multiply"].resume();
  block_diagonal_matrix_multiply(dX, Y, dY);
  timers[timerName + ".dY.multiply"].stop();
  timers[timerName + ".dY.subtract"].resume();
  dY -= R;
  timers[timerName + ".dY.subtract"].stop();
  timers[timerName + ".dY.cholesky"].resume();
  block_matrix_solve_with_cholesky(X_cholesky, dY);
  timers[timerName + ".dY.cholesky"].stop();
  timers[timerName + ".dY.symm"].resume();
  dY.symmetrize();
  dY *= Real(-1);
  dY *= El::BigFloat(-1);
  timers[timerName + ".dY.symm"].stop();
}
