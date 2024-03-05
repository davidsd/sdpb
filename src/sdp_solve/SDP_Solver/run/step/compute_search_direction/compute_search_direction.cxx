#include "sdp_solve/SDP_Solver/run/constraint_matrix_weighted_sum.hxx"

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

// C := alpha*A*B + beta*C
void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C);

// C := A B
inline void multiply(const Block_Diagonal_Matrix &A,
                     const Block_Diagonal_Matrix &B, Block_Diagonal_Matrix &C)
{
  scale_multiply_add(El::BigFloat(1), A, B, El::BigFloat(0), C);
}

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void cholesky_solve(const Block_Diagonal_Matrix &ACholesky,
                    Block_Diagonal_Matrix &X);

void compute_schur_RHS(const Block_Info &block_info, const SDP &sdp,
                       const Block_Vector &dual_residues,
                       const Block_Diagonal_Matrix &Z,
                       Block_Vector &dx);

void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);

void compute_search_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const bool &is_corrector_phase,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx,
  Block_Diagonal_Matrix &dX, Block_Vector &dy, Block_Diagonal_Matrix &dY)
{
  // R = beta mu I - X Y (predictor phase)
  // R = beta mu I - X Y - dX dY (corrector phase)
  Block_Diagonal_Matrix R(solver.X);

  scale_multiply_add(El::BigFloat(-1), solver.X, solver.Y, El::BigFloat(0), R);
  if(is_corrector_phase)
    {
      scale_multiply_add(El::BigFloat(-1), dX, dY, El::BigFloat(1), R);
    }
  R.add_diagonal(beta * mu);

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
  Block_Diagonal_Matrix Z(solver.X);
  multiply(solver.primal_residues, solver.Y, Z);
  Z -= R;
  cholesky_solve(X_cholesky, Z);
  Z.symmetrize();

  // dx[p] = -dual_residues[p] - Tr(A_p Z)
  // dy[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
  compute_schur_RHS(block_info, sdp, solver.dual_residues, Z, dx);
  dy=primal_residue_p;

  // Solve for dx, dy in-place
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal, Q, dx, dy);

  // dX = PrimalResidues + \sum_p A_p dx[p]
  constraint_matrix_weighted_sum(block_info, sdp, dx, dX);
  dX += solver.primal_residues;

  // dY = Symmetrize(X^{-1} (R - dX Y))
  multiply(dX, solver.Y, dY);
  dY -= R;
  cholesky_solve(X_cholesky, dY);
  dY.symmetrize();
  dY *= El::BigFloat(-1);
}
