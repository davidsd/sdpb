#include "../../lower_triangular_solve.hxx"
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
//

void compute_schur_complement(
  const SDP &sdp, const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Diagonal_Matrix &schur_complement);

void SDP_Solver::initialize_schur_complement_solver(
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  const std::vector<size_t> &block_dims,
  Block_Diagonal_Matrix &schur_complement_cholesky, Matrix &schur_off_diagonal,
  Block_Matrix &schur_off_diagonal_block, Matrix &Q,
  El::DistMatrix<El::BigFloat> &Q_elemental)
{
  // The Schur complement matrix S: a Block_Diagonal_Matrix with one
  // block for each 0 <= j < J.  SchurComplement.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  Block_Diagonal_Matrix schur_complement(block_dims);

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
  schur_off_diagonal = sdp.free_var_matrix;
  schur_off_diagonal_block = sdp.free_var_matrix_elemental;
  timers["initializeSchurComplementSolver.blockMatrixLowerTriangularSolve"]
    .resume();
  lower_triangular_solve(schur_complement_cholesky, schur_off_diagonal);
  lower_triangular_solve(schur_complement_cholesky, schur_off_diagonal_block);
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

  const size_t Q_size(schur_off_diagonal_block.blocks.at(0).Width());
  Zeros(Q_elemental, Q_size, Q_size);

  size_t schur_off_diagonal_total_rows(0);
  for(auto &block : schur_off_diagonal_block.blocks)
    {
      schur_off_diagonal_total_rows += block.Height();
    }
  El::DistMatrix<El::BigFloat> schur_off_diagonal_dist(
    schur_off_diagonal_total_rows, Q_size);

  size_t row_offset(0);
  for(auto &block : schur_off_diagonal_block.blocks)
    {
      El::DistMatrix<El::BigFloat> sub_matrix(
        El::View(schur_off_diagonal_dist, row_offset, 0, block.Height(),
                 block.Width()));
      El::Copy(block, sub_matrix);
      row_offset += block.Height();
    }

  // Here, SchurOffDiagonal = L'^{-1} B.
  //
  // UpperLeft(Q) = SchurOffDiagonal^T SchurOffDiagonal
  matrix_square_into_block(schur_off_diagonal, Q, 0, 0);

  // UpperRight(Q) = LowerLeft(Q)^T
  for(size_t c = 0; c < schur_off_diagonal.cols; c++)
    for(size_t r = schur_off_diagonal.cols; r < Q.rows; r++)
      {
        Q.elt(c, r) = Q.elt(r, c);
      }
  
  El::Syrk(El::UpperOrLowerNS::UPPER, El::OrientationNS::TRANSPOSE,
           El::BigFloat(1), schur_off_diagonal_dist, El::BigFloat(1),
           Q_elemental);
  El::MakeSymmetric(El::UpperOrLower::UPPER, Q_elemental);
  timers["initializeSchurComplementSolver.Qcomputation"].stop();

  timers["initializeSchurComplementSolver.LUDecomposition"].resume();
  timers["LUDecomposition.actualLU"].resume();

  Matrix Q_temp(Q);
  cholesky_decomposition(Q_temp, Q);
  Cholesky(El::UpperOrLowerNS::LOWER, Q_elemental);
  timers["LUDecomposition.actualLU"].stop();
  timers["initializeSchurComplementSolver.LUDecomposition"].stop();
}
