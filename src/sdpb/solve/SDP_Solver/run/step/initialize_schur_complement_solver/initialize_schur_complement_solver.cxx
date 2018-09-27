#include "../../../../SDP.hxx"
#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../Timers.hxx"

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
  Block_Diagonal_Matrix &schur_complement, Timers &timers);

void synchronize(const El::DistMatrix<El::BigFloat> &Q_group,
                 El::DistMatrix<El::BigFloat> &Q, Timers &timers);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y, const El::Grid &group_grid,
  Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q,
  Timers &timers)
{
  auto &initialize_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver"));
  // The Schur complement matrix S: a Block_Diagonal_Matrix with one
  // block for each 0 <= j < J.  SchurComplement.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  Block_Diagonal_Matrix schur_complement(
    block_info.schur_block_sizes, block_info.block_indices,
    block_info.schur_block_sizes.size(), group_grid);

  compute_schur_complement(block_info, bilinear_pairings_X_inv,
                           bilinear_pairings_Y, schur_complement, timers);

  auto &Q_computation_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.Qcomputation"));

  schur_off_diagonal.blocks.clear();
  schur_off_diagonal.blocks.reserve(schur_complement_cholesky.blocks.size());

  El::DistMatrix<El::BigFloat> Q_group(Q.Height(), Q.Width(), group_grid);
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

  synchronize(Q_group, Q, timers);
  Q_computation_timer.stop();

  auto &Cholesky_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver."
                         "Cholesky"));
  Cholesky(El::UpperOrLowerNS::UPPER, Q);
  Cholesky_timer.stop();
  initialize_timer.stop();
}
