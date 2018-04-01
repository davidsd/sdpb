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

void compute_schur_RHS(const SDP &sdp, const Vector &dualResidues,
                       const Block_Diagonal_Matrix &Z, const Vector &x,
                       Vector &r_x, Vector &r_y);

void SDP_Solver::compute_search_direction(const Real &beta, const Real &mu,
                                          const bool correctorPhase)
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

  timers[timerName + ".R.XY"].resume();
  block_diagonal_matrix_scale_multiply_add(-1, X, Y, 0, R);
  timers[timerName + ".R.XY"].stop();
  if(correctorPhase)
    {
      timers[timerName + ".R.dXdY"].resume();
      block_diagonal_matrix_scale_multiply_add(-1, dX, dY, 1, R);
      timers[timerName + ".R.dXdY"].stop();
    }

  timers[timerName + ".R.add"].resume();
  R.add_diagonal(beta * mu);
  timers[timerName + ".R.add"].stop();

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
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

  // r_x[p] = -dualResidues[p] - Tr(A_p Z)
  // r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
  // Here, dx = r_x, dy = r_y.
  timers[timerName + ".computeSchurRHS"].resume();
  compute_schur_RHS(sdp, dual_residues, Z, x, dx, dy);
  timers[timerName + ".computeSchurRHS"].stop();

  // Solve for dx, dy in-place
  timers[timerName + ".dxdy"].resume();
  solve_schur_complement_equation(dx, dy);
  timers[timerName + ".dxdy"].stop();

  // dX = PrimalResidues + \sum_p A_p dx[p]
  timers[timerName + ".dX.weightedSum"].resume();
  constraint_matrix_weighted_sum(sdp, dx, dX);
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
  dY *= -1;
  timers[timerName + ".dY.symm"].stop();
}
