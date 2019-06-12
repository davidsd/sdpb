#include "../../../../SDP.hxx"
#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../../Timers.hxx"

void initialize_Q_group(const SDP &sdp, const Block_Info &block_info,
                        const Block_Diagonal_Matrix &S, Block_Matrix &L_inv_B,
                        Block_Diagonal_Matrix &L,
                        El::DistMatrix<El::BigFloat> &Q_group,
                        Block_Diagonal_Matrix &eigenvectors, Timers &timers)
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

  L_inv_B.blocks.clear();
  L_inv_B.blocks.reserve(L.blocks.size());
  eigenvectors.blocks.clear();
  eigenvectors.blocks.reserve(S.blocks.size());

  El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
  /// The default number of iterations is 40.  That is sometimes
  /// not enough, so we bump it up significantly.
  hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 512;

  for(size_t block = 0; block < L.blocks.size(); block++)
    {
      auto &cholesky_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.cholesky_"
        + std::to_string(block_info.block_indices[block])));
      L.blocks[block] = S.blocks[block];

      El::DistMatrix<El::BigFloat> eigenvalues(S.blocks[block].Grid());
      eigenvectors.blocks.emplace_back(S.blocks[block].Grid());
      auto block_copy(S.blocks[block]);
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff
        = S.blocks[block].Height() / 2 + 1;
      El::HermitianEig(El::UpperOrLowerNS::LOWER, block_copy, eigenvalues,
                       eigenvectors.blocks.back(), hermitian_eig_ctrl);

      El::BigFloat max_eigenvalue(El::Max(eigenvalues));
      const El::BigFloat eigenvalue_cutoff(
        1e-50 * max_eigenvalue); // FIXME: This is totally arbitrary
      for(int64_t row = 0; row < eigenvalues.LocalHeight(); ++row)
        for(int64_t column = 0; column < eigenvalues.LocalWidth(); ++column)
          {
            El::BigFloat local_value(eigenvalues.GetLocal(row, column));
            if(local_value < eigenvalue_cutoff)
              {
                eigenvalues.SetLocal(row, column, El::BigFloat(0));
              }
            else
              {
                eigenvalues.SetLocal(row, column, El::Sqrt(local_value));
              }
          }

      Cholesky(El::UpperOrLowerNS::LOWER, L.blocks[block]);
      cholesky_timer.stop();

      // L_inv_B = L^-1 B
      auto &solve_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.solve_"
        + std::to_string(block_info.block_indices[block])));

      L_inv_B.blocks.push_back(sdp.B.blocks[block]);
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), L.blocks[block], L_inv_B.blocks[block]);

      solve_timer.stop();

      // Q = L_inv_B^T L_inv_B
      auto &syrk_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.syrk_"
        + std::to_string(block_info.block_indices[block])));
      El::DistMatrix<El::BigFloat> Q_group_view(
        El::View(Q_group, 0, 0, L_inv_B.blocks[block].Width(),
                 L_inv_B.blocks[block].Width()));
      El::Syrk(El::UpperOrLowerNS::UPPER, El::OrientationNS::TRANSPOSE,
               El::BigFloat(1), L_inv_B.blocks[block], El::BigFloat(1),
               Q_group_view);
      syrk_timer.stop();
    }
}
