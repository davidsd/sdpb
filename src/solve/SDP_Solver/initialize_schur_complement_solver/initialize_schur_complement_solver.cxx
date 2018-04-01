#include "../../SDP_Solver.hxx"
#include "../../../Timers.hxx"

// Compute the quantities needed to solve the Schur complement
// equation
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {r, s}
//
// (where S = SchurComplement, B = FreeVarMatrix), using the method
// described in the manual:
//
// - Compute S using BilinearPairingsXInv and BilinearPairingsY.
//
// - Compute the Cholesky decomposition S' = L' L'^T.
//
// - Form B' = (B U) and compute
//
//   - SchurOffDiagonal = L'^{-1} B
//   - L'^{-1} U
//   - Q = (L'^{-1} B')^T (L'^{-1} B') - {{0, 0}, {0, 1}}
//
// - Compute the LU decomposition of Q.
//
// This data is sufficient to efficiently solve the above equation for
// a given r,s.
//
// Inputs:
// - BilinearPairingsXInv, BilinearPairingsY (these are members of
//   SDPSolver, but we include them as arguments to emphasize that
//   they must be computed first)
// Workspace (members of SDPSolver which are modified by this method
// and not used later):
// - SchurComplement
// Outputs (members of SDPSolver which are modified by this method and
// used later):
// - SchurComplementCholesky
// - SchurOffDiagonal
// - Q, Qpivots
//

void compute_schur_complement(
  const SDP &sdp, const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Diagonal_Matrix &schur_complement);

void SDP_Solver::initialize_schur_complement_solver(
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y)
{
  compute_schur_complement(sdp, bilinear_pairings_X_inv, bilinear_pairings_Y,
                           schur_complement);

  // compute SchurComplementCholesky = L', where
  //
  //   L' L'^T = S'
  //
  timers["initializeSchurComplementSolver.choleskyDecomposition"].resume();
  cholesky_decomposition(schur_complement, schur_complement_cholesky);
  timers["initializeSchurComplementSolver.choleskyDecomposition"].stop();

  // SchurOffDiagonal = L'^{-1} FreeVarMatrix
  schur_off_diagonal.copy_from(sdp.free_var_matrix);
  timers["initializeSchurComplementSolver.blockMatrixLowerTriangularSolve"]
    .resume();
  block_matrix_lower_triangular_solve(schur_complement_cholesky,
                                      schur_off_diagonal);
  timers["initializeSchurComplementSolver.blockMatrixLowerTriangularSolve"]
    .stop();

  // total number of columns in the off-diagonal part B' = (B U)
  // (currently just B; will accumulate the rest shortly)
  int off_diagonal_columns = schur_off_diagonal.cols;

  // Next, we compute
  //
  //   Q = (L'^{-1} B')^T (L'^{-1} B') - {{0, 0}, {0, 1}}
  //
  // Where B' = (B U).  We think of Q as containing four blocks called
  // Upper/Lower-Left/Right.

  timers["initializeSchurComplementSolver.Qcomputation"].resume();
  // Set the dimensions of Q
  Q.resize(off_diagonal_columns, off_diagonal_columns);
  Q.set_zero();

  // Here, SchurOffDiagonal = L'^{-1} B.
  //
  // UpperLeft(Q) = SchurOffDiagonal^T SchurOffDiagonal
  matrix_square_into_block(schur_off_diagonal, Q, 0, 0);

  // UpperRight(Q) = LowerLeft(Q)^T
  for(int c = 0; c < schur_off_diagonal.cols; c++)
    for(int r = schur_off_diagonal.cols; r < Q.rows; r++)
      Q.elt(c, r) = Q.elt(r, c);
  timers["initializeSchurComplementSolver.Qcomputation"].stop();

  // Resize Qpivots appropriately and compute the LU decomposition of Q

  timers["initializeSchurComplementSolver.LUDecomposition"].resume();
  timers["LUDecomposition.resizing"].resume();
  Q_pivots.resize(Q.rows);
  timers["LUDecomposition.resizing"].stop();
  timers["LUDecomposition.actualLU"].resume();
  LU_decomposition(Q, Q_pivots);
  timers["LUDecomposition.actualLU"].stop();
  timers["initializeSchurComplementSolver.LUDecomposition"].stop();
}
