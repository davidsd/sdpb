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

void compute_schur_RHS(const SDP &sdp,
                       const Block_Vector &dual_residues_elemental,
                       const Block_Diagonal_Matrix &Z,
                       const Block_Vector &x_elemental,
                       Block_Vector &r_x_elemental,
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
  const El::DistMatrix<El::BigFloat> &Q_elemental, Block_Vector &dx_elemental,
  Block_Diagonal_Matrix &dX, El::DistMatrix<El::BigFloat> &dy_elemental,
  Block_Diagonal_Matrix &dY)
{
  // R = beta mu I - X Y (predictor phase)
  // R = beta mu I - X Y - dX dY (corrector phase)
  Block_Diagonal_Matrix R(X);

  block_diagonal_matrix_scale_multiply_add(El::BigFloat(-1), X, Y,
                                           El::BigFloat(0), R);
  if(correctorPhase)
    {
      block_diagonal_matrix_scale_multiply_add(El::BigFloat(-1), dX, dY,
                                               El::BigFloat(1), R);
    }

  R.add_diagonal(beta_elemental * mu_elemental);

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))

  // FIXME: It has to have the same shape as X, but does not need to
  // be a copy.
  Block_Diagonal_Matrix Z(X);

  block_diagonal_matrix_multiply(primal_residues, Y, Z);
  Z -= R;
  block_matrix_solve_with_cholesky(X_cholesky, Z);
  Z.symmetrize();

  // r_x[p] = -dual_residues[p] - Tr(A_p Z)
  // r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
  // Here, dx = r_x, dy = r_y.
  compute_schur_RHS(sdp, dual_residues_elemental, Z, x_elemental, dx_elemental,
                    dy_elemental);

  // Solve for dx, dy in-place
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal_elemental, Q_elemental,
                                  dx_elemental, dy_elemental);

  // dX = PrimalResidues + \sum_p A_p dx[p]
  constraint_matrix_weighted_sum(sdp, dx_elemental, dX);
  dX += primal_residues;

  // dY = Symmetrize(X^{-1} (R - dX Y))
  block_diagonal_matrix_multiply(dX, Y, dY);
  dY -= R;
  block_matrix_solve_with_cholesky(X_cholesky, dY);
  dY.symmetrize();
  dY *= El::BigFloat(-1);
}
