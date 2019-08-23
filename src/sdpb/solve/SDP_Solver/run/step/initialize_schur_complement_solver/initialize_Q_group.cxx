#include "../../../../SDP.hxx"
#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../../Timers.hxx"

void initialize_Q_group(const SDP &sdp, const Block_Info &block_info,
                        const Block_Diagonal_Matrix &schur_complement,
                        Block_Matrix &schur_off_diagonal,
                        Block_Diagonal_Matrix &schur_complement_cholesky,
                        El::DistMatrix<El::BigFloat> &Q_group, Timers &timers)
{
  // Explicitly deallocate the lower half of Q_group.  This
  // significantly reduces the total amount of memory required.
  El::Matrix<El::BigFloat> &local(Q_group.Matrix());
  for(int64_t row = 0; row < Q_group.Height(); ++row)
    for(int64_t column = 0; column < row; ++column)
      {
        if(Q_group.IsLocal(row, column))
          {
            mpf_clear(local(Q_group.LocalRow(row), Q_group.LocalCol(column))
                        .gmp_float.get_mpf_t());
            local(Q_group.LocalRow(row), Q_group.LocalCol(column))
              .gmp_float.get_mpf_t()[0]
              ._mp_d
              = nullptr;
          }
      }

  schur_off_diagonal.blocks.clear();
  schur_off_diagonal.blocks.reserve(schur_complement_cholesky.blocks.size());

  El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
  /// The default number of iterations is 40.  That is sometimes
  /// not enough, so we bump it up significantly.
  hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 512;
  
  for(size_t block = 0; block < schur_complement_cholesky.blocks.size();
      block++)
    {
      auto &cholesky_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.cholesky_"
        + std::to_string(block_info.block_indices[block])));

      {
      El::DistMatrix<El::BigFloat> eigenvalues(schur_complement.blocks[block].Grid());
      auto block_copy(schur_complement.blocks[block]);
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff
        = schur_complement.blocks[block].Height() / 2 + 1;
      El::HermitianEig(El::UpperOrLowerNS::LOWER, block_copy, eigenvalues,
                       hermitian_eig_ctrl);
      std::ofstream eigvals("eigvals_" + std::to_string(block_info.block_indices[block]) + ".txt");
      El::Print(eigenvalues,"","\n",eigvals);
      // std::cout << "eigvals: " << El::mpi::Rank() << " "
      //           << El::Max(eigenvalues)/El::Min(eigenvalues) << " "
      //           << "\n";
      }

      
      schur_complement_cholesky.blocks[block] = schur_complement.blocks[block];
      Cholesky(El::UpperOrLowerNS::LOWER,
               schur_complement_cholesky.blocks[block]);
      cholesky_timer.stop();

      // SchurOffDiagonal = L'^{-1} FreeVarMatrix
      auto &solve_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.solve_"
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
        "run.step.initializeSchurComplementSolver.Q.syrk_"
        + std::to_string(block_info.block_indices[block])));
      El::DistMatrix<El::BigFloat> Q_group_view(
        El::View(Q_group, 0, 0, schur_off_diagonal.blocks[block].Width(),
                 schur_off_diagonal.blocks[block].Width()));
      El::Syrk(El::UpperOrLowerNS::UPPER, El::OrientationNS::TRANSPOSE,
               El::BigFloat(1), schur_off_diagonal.blocks[block],
               El::BigFloat(1), Q_group_view);
      syrk_timer.stop();
    }
}
