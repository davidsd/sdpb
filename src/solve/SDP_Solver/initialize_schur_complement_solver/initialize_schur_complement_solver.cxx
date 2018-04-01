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
// - Stabilize S by adding a low-rank update S' = S + U U^T and
//   compute the Cholesky decomposition S' = L' L'^T.
//
// - Form B' = (B U) and compute
//
//   - SchurOffDiagonal = L'^{-1} B
//   - L'^{-1} U (this is stored implicitly in the stabilizeBlock* variables)
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
// - choleskyStabilizeThreshold: the real constant \theta used in
//   stabilizing S
// Workspace (members of SDPSolver which are modified by this method
// and not used later):
// - SchurComplement
// - schurStabilizeIndices
// - schurStabilizeLambdas
// Outputs (members of SDPSolver which are modified by this method and
// used later):
// - SchurComplementCholesky
// - SchurOffDiagonal
// - stabilizeBlockIndices
// - stabilizeBlockUpdateRow
// - stabilizeBlockUpdateColumn
// - stabilizeBlocks
// - Q, Qpivots
//

void compute_schur_complement(
  const SDP &sdp, const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Diagonal_Matrix &schur_complement);

void SDP_Solver::initialize_schur_complement_solver(
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  const Real &cholesky_stabilize_threshold)
{
  compute_schur_complement(sdp, bilinear_pairings_X_inv, bilinear_pairings_Y,
                           schur_complement);

  // compute SchurComplementCholesky = L', where
  //
  //   L' L'^T = S' = SchurComplement + U U^T
  //
  // Here, the 'update' matrix U has columns given by
  //
  //   U = ( Lambda_{p_1} e_{p_1}, ..., Lambda_{p_M} e_{p_M} )
  //
  // where e_p is a unit vector in the p-th direction and the
  // Lambda_{p_m} are constants. The p_i are given by
  // schurStabilizeIndices and the corresponding Lambda_i by
  // schurStabilizeLambdas.
  //
  timers["initializeSchurComplementSolver.choleskyDecompositionStabilized"]
    .resume();
  cholesky_decomposition_stabilized(
    schur_complement, schur_complement_cholesky, schur_stabilize_indices,
    schur_stabilize_lambdas, cast2double(cholesky_stabilize_threshold));
  timers["initializeSchurComplementSolver.choleskyDecompositionStabilized"]
    .stop();

  // SchurOffDiagonal = L'^{-1} FreeVarMatrix
  schur_off_diagonal.copy_from(sdp.free_var_matrix);
  timers["initializeSchurComplementSolver.blockMatrixLowerTriangularSolve"]
    .resume();
  block_matrix_lower_triangular_solve(schur_complement_cholesky,
                                      schur_off_diagonal);
  timers["initializeSchurComplementSolver.blockMatrixLowerTriangularSolve"]
    .stop();

  // Next we compute L'^{-1} U, which is stored implicitly in terms of
  // its nonzero submatrices in the stabilizeBlock* variables.

  // total number of columns in the off-diagonal part B' = (B U)
  // (currently just B; will accumulate the rest shortly)
  int off_diagonal_columns = schur_off_diagonal.cols;

  stabilize_block_indices.resize(0);
  stabilize_block_update_Row.resize(0);
  stabilize_block_update_column.resize(0);

  // j runs over blocks of SchurComplement
  timers["initializeSchurComplementSolver.stabilizeMagic"].resume();
  for(size_t j = 0; j < schur_complement.blocks.size(); j++)
    {
      if(schur_stabilize_indices[j].size() > 0)
        {
          // the j-th block of S contains stabilized directions. We have a
          // block submatrix
          //
          //   U_j = stabilizeBlocks[j]
          //
          // of U for each such j.
          stabilize_block_indices.push_back(j);

          // index of the first stabilized direction within the j-th block
          int start_index = schur_stabilize_indices[j][0];

          // set the dimensions of U_j
          stabilize_blocks[j].resize(schur_complement.blocks[j].rows
                                       - start_index,
                                     schur_stabilize_indices[j].size());
          // set U_j = 0
          stabilize_blocks[j].set_zero();
          // for each column of U_j add Lambda_p in the row (p - startIndex)
          for(size_t c = 0; c < schur_stabilize_indices[j].size(); c++)
            {
              int r = schur_stabilize_indices[j][c] - start_index;
              stabilize_blocks[j].elt(r, c) = schur_stabilize_lambdas[j][c];
            }

          // append the row of U corresponding to the top-left of U_j
          stabilize_block_update_Row.push_back(
            schur_complement.block_start_indices[j] + start_index);
          // append the column of U corresponding to the top-left of U_j
          stabilize_block_update_column.push_back(off_diagonal_columns);
          // update the number of off-diagonal columns
          off_diagonal_columns += stabilize_blocks[j].cols;
        }
    }
  timers["initializeSchurComplementSolver.stabilizeMagic"].stop();

  // Set U = L'^{-1} U
  //
  // We do this by modifying the blocks U_j = stabilizeBlocks[j]
  // in-place, multiplying by the inverse of the appropriate submatrix
  // of SchurComplementCholesky.  We henceforth refer to L'^{-1} U as
  // V to avoid confusion.
  timers["initializeSchurComplementSolver.LinvU"].resume();
  for(unsigned int j = 0; j < stabilize_block_indices.size(); j++)
    {
      int b = stabilize_block_indices[j];
      int startIndex = schur_stabilize_indices[b][0];
      Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal",
            stabilize_blocks[b].rows, stabilize_blocks[b].cols, 1,
            &schur_complement_cholesky.blocks[b].elt(startIndex, startIndex),
            schur_complement_cholesky.blocks[b].rows,
            &stabilize_blocks[b].elt(0, 0), stabilize_blocks[b].rows);
    }
  timers["initializeSchurComplementSolver.LinvU"].stop();

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

  // Here, stabilizeBlocks contains the blocks of V = L'^{-1} U.
  //
  // LowerRight(Q) = V^T V - 1
  for(unsigned int j = 0; j < stabilize_block_indices.size(); j++)
    {
      int b = stabilize_block_indices[j];
      int c = stabilize_block_update_column[j];
      matrix_square_into_block(stabilize_blocks[b], Q, c, c);

      // subtract the identity matrix from this block
      for(int i = c; i < c + stabilize_blocks[b].cols; i++)
        Q.elt(i, i) -= 1;
    }

  // LowerLeft(Q) = V^T SchurOffDiagonal
  timers["Qcomputation.nonGPU"].resume();
  for(unsigned int j = 0; j < stabilize_block_indices.size(); j++)
    {
      int b = stabilize_block_indices[j];
      int p = stabilize_block_update_Row[j];
      int r = stabilize_block_update_column[j];

      Rgemm("Transpose", "NoTranspose", stabilize_blocks[b].cols,
            schur_off_diagonal.cols, stabilize_blocks[b].rows, 1,
            &stabilize_blocks[b].elements[0], stabilize_blocks[b].rows,
            &schur_off_diagonal.elt(p, 0), schur_off_diagonal.rows, 0,
            &Q.elt(r, 0), Q.rows);
    }
  timers["Qcomputation.nonGPU"].stop();

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
