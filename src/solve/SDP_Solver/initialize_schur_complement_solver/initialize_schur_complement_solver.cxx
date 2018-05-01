#include "../../SDP_Solver.hxx"
#include "../../../Timers.hxx"

#include <El.hpp>

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
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  const std::vector<size_t> &block_dims,
  Block_Diagonal_Matrix &schur_complement_cholesky, Matrix &schur_off_diagonal,
  Matrix &Q)
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

  // const El::Grid grid(El::mpi::COMM_WORLD);
  // El::DistMatrix<El::BigFloat> schur_off_diagonal_elemental(grid),
  //   Q_elemental(grid);
  // schur_off_diagonal_elemental.Resize(schur_off_diagonal.rows,
  //                                     schur_off_diagonal.cols);
  // Zeros(Q_elemental, schur_off_diagonal.cols, schur_off_diagonal.cols);

  // const El::Int local_height(schur_off_diagonal_elemental.LocalHeight()),
  //   local_width(schur_off_diagonal_elemental.LocalWidth());

  // for(size_t row = 0; row < local_height; ++row)
  //   {
  //     El::Int global_row(schur_off_diagonal_elemental.GlobalRow(row));
  //     for(size_t column = 0; column < local_width; ++column)
  //       {
  //         El::Int global_column(
  //           schur_off_diagonal_elemental.GlobalCol(column));

  //         std::stringstream ss;
  //         ss.precision(1024 * 0.31 + 5);
  //         ss << schur_off_diagonal.elt(global_row, global_column);
  //         schur_off_diagonal_elemental.SetLocal(row, column,
  //                                               El::BigFloat(ss.str(), 10));
  //       }
  //   }
  // El::Syrk(El::UPPER, El::TRANSPOSE, El::BigFloat(1),
  //          schur_off_diagonal_elemental, El::BigFloat(1), Q_elemental);

  // const El::Int Q_height(Q_elemental.LocalHeight()),
  //   Q_width(Q_elemental.LocalWidth());

  // for(size_t row = 0; row < Q_height; ++row)
  //   {
  //     El::Int global_row(Q_elemental.GlobalRow(row));
  //     for(size_t column = 0; column < Q_width; ++column)
  //       {
  //         El::Int global_column(Q_elemental.GlobalCol(column));

  //         std::stringstream ss;
  //         ss.precision(1024 * 0.31 + 5);
  //         ss << Q_elemental.GetLocal(row, column);
  //         Q.elt(global_row, global_column) = ss.str();
  //       }
  //   }

  // std::cout.precision(1024 * 0.31 + 5);
  // std::cout << "Qe: "
  //           // << local_height << " " << local_width << " "
  //           << Q_height << " " << Q_width << " " << Q.rows << " " << Q.cols
  //           << " " << Q_elemental.GetLocal(0, 0) << " "
  //           << "\n";
  // std::cout << "Qe: "
  //           // << local_height << " " << local_width << " "
  //           << Q_height << " " << Q_width << " " << Q.rows << " " << Q.cols
  //           << " " << Q.elt(0, 0) << " "
  //           << "\n";
  // std::cout << El::mpfr::Precision() << " "
  //           << schur_off_diagonal.elt(0,0).get_prec() << " "
  //           << Q.elt(0, 0).get_prec() << " "
  //            << "\n";

  // // std::cout << "Qm: " << Q_height << " " << Q_width << " " << Q.rows << "
  // "
  // //           << Q.cols << " " << Q.elt(0, 0) << " "
  // //           << "\n";

  // UpperRight(Q) = LowerLeft(Q)^T
  for(size_t c = 0; c < schur_off_diagonal.cols; c++)
    for(size_t r = schur_off_diagonal.cols; r < Q.rows; r++)
      {
        Q.elt(c, r) = Q.elt(r, c);
      }
  timers["initializeSchurComplementSolver.Qcomputation"].stop();

  // Resize Qpivots appropriately and compute the LU decomposition of Q

  timers["initializeSchurComplementSolver.LUDecomposition"].resume();
  timers["LUDecomposition.actualLU"].resume();

  Matrix Q_temp(Q);
  cholesky_decomposition(Q_temp,Q);
  timers["LUDecomposition.actualLU"].stop();
  timers["initializeSchurComplementSolver.LUDecomposition"].stop();
}
