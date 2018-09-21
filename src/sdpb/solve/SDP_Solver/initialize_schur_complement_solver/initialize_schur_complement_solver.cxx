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
  const Block_Info &block_info,
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Diagonal_Matrix &schur_complement);

void SDP_Solver::initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y, const El::Grid &block_grid,
  const bool &debug, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q,
  Timers &timers)
{
  auto &schur_complement_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.schur_complement"));
  // The Schur complement matrix S: a Block_Diagonal_Matrix with one
  // block for each 0 <= j < J.  SchurComplement.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  Block_Diagonal_Matrix schur_complement(
    block_info.schur_block_sizes, block_info.block_indices,
    block_info.schur_block_sizes.size(), block_grid);

  compute_schur_complement(block_info, bilinear_pairings_X_inv,
                           bilinear_pairings_Y, schur_complement);
  schur_complement_timer.stop();

  if(debug)
    {
      El::Output(El::mpi::Rank(), " run.step.initializeSchurComplementSolver."
                                  "Qcomputation");
    }
  auto &Q_computation_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.Qcomputation"));

  schur_off_diagonal.blocks.clear();
  schur_off_diagonal.blocks.reserve(schur_complement_cholesky.blocks.size());

  El::DistMatrix<El::BigFloat> Q_group(Q.Height(), Q.Width(), block_grid);
  El::Zero(Q_group);
  for(size_t block = 0; block < schur_complement_cholesky.blocks.size();
      block++)
    {
      auto &cholesky_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Qcomputation.cholesky_"
        + std::to_string(block_info.block_indices[block])));
      schur_complement_cholesky.blocks[block] = schur_complement.blocks[block];
      Cholesky(El::UpperOrLowerNS::LOWER,
               schur_complement_cholesky.blocks[block]);
      cholesky_timer.stop();

      // SchurOffDiagonal = L'^{-1} FreeVarMatrix
      auto &solve_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver."
        "Qcomputation.solve_"
        + std::to_string(block_info.block_indices[block])));

      schur_off_diagonal.blocks.push_back(sdp.free_var_matrix.blocks[block]);
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), schur_complement_cholesky.blocks[block],
               schur_off_diagonal.blocks[block]);

      solve_timer.stop();

      // Next, we compute
      //
      //   Q = (L'^{-1} B')^T (L'^{-1} B') - {{0, 0}, {0, 1}}
      //
      // Where B' = (B U).  We think of Q as containing four blocks
      // called Upper/Lower-Left/Right.

      auto &syrk_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver."
        "Qcomputation.syrk_"
        + std::to_string(block_info.block_indices[block])));
      El::DistMatrix<El::BigFloat> Q_group_view(
        El::View(Q_group, 0, 0, schur_off_diagonal.blocks[block].Width(),
                 schur_off_diagonal.blocks[block].Width()));
      El::Syrk(El::UpperOrLowerNS::UPPER, El::OrientationNS::TRANSPOSE,
               El::BigFloat(1), schur_off_diagonal.blocks[block],
               El::BigFloat(1), Q_group_view);
      syrk_timer.stop();
    }

  // Synchronize the results back to the global Q.

  auto &synchronize_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.Qcomputation."
    "synchronize"));
  // Optimize for when Q_group is on a single processor
  if(Q_group.Grid().Size() == 1)
    {
      El::AllReduce(Q_group, El::mpi::COMM_WORLD);
      El::MakeSymmetric(El::UpperOrLower::UPPER, Q_group);

      for(int64_t row = 0; row < Q.LocalHeight(); ++row)
        {
          int64_t global_row(Q.GlobalRow(row));
          for(int64_t column = 0; column < Q.LocalWidth(); ++column)
            {
              int64_t global_column(Q.GlobalCol(column));
              Q.SetLocal(row, column,
                         Q_group.GetLocal(global_row, global_column));
            }
        }
    }
  else
    {
      El::Matrix<El::BigFloat> Q_local(Q.Height(), Q.Width());
      El::Zero(Q_local);

      for(int64_t row = 0; row < Q_group.LocalHeight(); ++row)
        {
          int64_t global_row(Q_group.GlobalRow(row));
          for(int64_t column = 0; column < Q_group.LocalWidth(); ++column)
            {
              int64_t global_column(Q_group.GlobalCol(column));
              Q_local(global_row, global_column)
                = Q_group.GetLocal(row, column);
            }
        }

      El::AllReduce(Q_local, El::mpi::COMM_WORLD);
      El::MakeSymmetric(El::UpperOrLower::UPPER, Q_local);

      for(int64_t row = 0; row < Q.LocalHeight(); ++row)
        {
          int64_t global_row(Q.GlobalRow(row));
          for(int64_t column = 0; column < Q.LocalWidth(); ++column)
            {
              int64_t global_column(Q.GlobalCol(column));
              Q.SetLocal(row, column, Q_local(global_row, global_column));
            }
        }
    }
  synchronize_timer.stop();
  Q_computation_timer.stop();

  if(debug)
    {
      El::Output(El::mpi::Rank(),
                 " run.step.initializeSchurComplementSolver.Cholesky");
    }
  auto &Cholesky_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver."
                         "Cholesky"));
  Cholesky(El::UpperOrLowerNS::LOWER, Q);
  Cholesky_timer.stop();
}
