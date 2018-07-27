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
  const std::vector<size_t> &block_sizes,
  const std::vector<size_t> &block_indices, const size_t &num_schur_blocks,
  const El::Grid &block_grid, const std::vector<size_t> &schur_offsets,
  const bool &debug, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal_block, El::DistMatrix<El::BigFloat> &Q)
{
  // The Schur complement matrix S: a Block_Diagonal_Matrix with one
  // block for each 0 <= j < J.  SchurComplement.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  Block_Diagonal_Matrix schur_complement(block_sizes, block_indices,
                                         num_schur_blocks, block_grid);

  compute_schur_complement(sdp, bilinear_pairings_X_inv, bilinear_pairings_Y,
                           schur_complement);

  // compute SchurComplementCholesky = L', where
  //
  //   L' L'^T = S'
  //
  if(debug)
    {
      El::Output(El::mpi::Rank(),
        " run.step.initializeSchurComplementSolver.choleskyDecomposition");
    }
  timers["run.step.initializeSchurComplementSolver.choleskyDecomposition"]
    .resume();
  cholesky_decomposition(schur_complement, schur_complement_cholesky);
  timers["run.step.initializeSchurComplementSolver.choleskyDecomposition"]
    .stop();

  // SchurOffDiagonal = L'^{-1} FreeVarMatrix
  schur_off_diagonal_block = sdp.free_var_matrix;
  if(debug)
    {
      El::Output(El::mpi::Rank()," run.step.initializeSchurComplementSolver."
                         "blockMatrixLowerTriangularSolve");
    }
  timers["run.step.initializeSchurComplementSolver."
         "blockMatrixLowerTriangularSolve"]
    .resume();
  lower_triangular_solve(schur_complement_cholesky, schur_off_diagonal_block);
  timers["run.step.initializeSchurComplementSolver."
         "blockMatrixLowerTriangularSolve"]
    .stop();

  // Next, we compute
  //
  //   Q = (L'^{-1} B')^T (L'^{-1} B') - {{0, 0}, {0, 1}}
  //
  // Where B' = (B U).  We think of Q as containing four blocks called
  // Upper/Lower-Left/Right.

  if(debug)
    {
      El::Output(El::mpi::Rank(),
        " run.step.initializeSchurComplementSolver.Qcomputation");
    }
  timers["run.step.initializeSchurComplementSolver.Qcomputation"].resume();
  El::Zero(Q);
  El::DistMatrix<El::BigFloat> schur_off_diagonal_dist(schur_offsets.back(),
                                                       Q.Width());
  Zero(schur_off_diagonal_dist);
  {
    size_t num_updates(0);
    for(auto &block : schur_off_diagonal_block.blocks)
      {
        num_updates += block.Height() * block.Width();
      }
    schur_off_diagonal_dist.Reserve(num_updates);
  }
  {
    auto schur_off_diagonal_block_block(
      schur_off_diagonal_block.blocks.begin());
    for(auto &block_index : block_indices)
      {
        for(int64_t row = 0; row < schur_off_diagonal_block_block->Height();
            ++row)
          for(int64_t column = 0;
              column < schur_off_diagonal_block_block->Width(); ++column)
            {
              schur_off_diagonal_dist.QueueUpdate(
                row + schur_offsets.at(block_index), column,
                schur_off_diagonal_block_block->GetLocal(row, column));
            }
        ++schur_off_diagonal_block_block;
      }
    schur_off_diagonal_dist.ProcessQueues();
  }

  // Here, SchurOffDiagonal = L'^{-1} B.
  //
  // UpperLeft(Q) = SchurOffDiagonal^T SchurOffDiagonal
  if(debug)
    {
      El::Output(El::mpi::Rank()," run.step.initializeSchurComplementSolver.Syrk");
    }
  El::Syrk(El::UpperOrLowerNS::UPPER, El::OrientationNS::TRANSPOSE,
           El::BigFloat(1), schur_off_diagonal_dist, El::BigFloat(1), Q);

  El::MakeSymmetric(El::UpperOrLower::UPPER, Q);
  timers["run.step.initializeSchurComplementSolver.Qcomputation"].stop();

  if(debug)
    {
      El::Output(El::mpi::Rank(),
        " run.step.initializeSchurComplementSolver.LUDecomposition");
    }
  timers["run.step.initializeSchurComplementSolver.LUDecomposition"].resume();
  Cholesky(El::UpperOrLowerNS::LOWER, Q);
  timers["run.step.initializeSchurComplementSolver.LUDecomposition"].stop();
}
