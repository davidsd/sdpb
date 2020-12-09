#include "../../../../SDP.hxx"
#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../Timers.hxx"

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

  // size_t num_rows(sdp.free_var_matrix.blocks.front().Width());
  // for(auto &block: sdp.free_var_matrix.blocks)
  //   {
  //     num_rows+=block.Height();
  //   }
  // // std::cout << "rows: " << num_rows << "\n";

  // El::Matrix<double> big_matrix(num_rows,num_rows);
  // El::Zero(big_matrix);
  // {
  //   int64_t row(0);
  //   for(auto &block: schur_complement.blocks)
  //     {
  //       // El::Print(block,"block");
  //       // std::cout << "\n";
  //       for(int64_t block_row(0); block_row!=block.Height(); ++block_row)
  //         {
  //           for(int64_t block_column(0); block_column!=block.Width(); ++block_column)
  //             {
  //               big_matrix(row+block_row, row+block_column) = double(block.Get(block_row,block_column));
  //             }
  //         }
  //       row+=block.Height();
  //     }
  // }
  // {
  //   int64_t row(0);
  //   const int64_t column_offset(num_rows - sdp.free_var_matrix.blocks.front().Width());
  //   for(auto &block: sdp.free_var_matrix.blocks)
  //     {
  //       for(int64_t block_row(0); block_row!=block.Height(); ++block_row)
  //         {
  //           for(int64_t block_column(0); block_column!=block.Width(); ++block_column)
  //             {
  //               big_matrix(row+block_row, column_offset+block_column) =
  //                 double(block.Get(block_row,block_column));
  //               big_matrix(column_offset+block_column, row+block_row) =
  //                 big_matrix(row+block_row, column_offset+block_column);
  //             }
  //         }
  //       row+=block.Height();
  //     }
  // }
  // {
  //   // El::Print(big_matrix,"big_matrix");
  //   // std::cout << "\n";
  //   // exit(0);
  //   /// There is a bug in El::HermitianEig when there is more than
  //   /// one level of recursion when computing eigenvalues.  One fix
  //   /// is to increase the cutoff so that there is no more than one
  //   /// level of recursion.

  //   /// An alternate workaround is to compute both eigenvalues and
  //   /// eigenvectors, but that seems to be significantly slower.
  //   El::HermitianEigCtrl<double> hermitian_eig_ctrl;
  //   hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff = num_rows / 2 + 1;

  //   /// The default number of iterations is 40.  That is sometimes
  //   /// not enough, so we bump it up significantly.
  //   hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 400;

  //   El::Matrix<double> eigenvalues;
  //   El::HermitianEig(El::UpperOrLowerNS::LOWER, big_matrix, eigenvalues,
  //                    hermitian_eig_ctrl);

  //   double max(0), min(std::numeric_limits<double>::max());

  //   for(int64_t row(0); row!=eigenvalues.Height(); ++row)
  //     {
  //       max=std::max(max,std::abs(eigenvalues(row,0)));
  //       min=std::min(min,std::abs(eigenvalues(row,0)));
  //     }
  //   std::cout << "T Condition: "
  //             << (max / min) << " "
  //             << max << " "
  //             << min << "\n";
  // }
  
  El::BigFloat S_max(0), S_min(std::numeric_limits<double>::max()),
    S_condition(0);

  for(size_t block = 0; block < schur_complement_cholesky.blocks.size();
      ++block)
    {
      auto &cholesky_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.cholesky_"
        + std::to_string(block_info.block_indices[block])));
      schur_complement_cholesky.blocks[block] = schur_complement.blocks[block];

      if(schur_complement.blocks[block].Height() > 1)
      {
        /// There is a bug in El::HermitianEig when there is more than
        /// one level of recursion when computing eigenvalues.  One fix
        /// is to increase the cutoff so that there is no more than one
        /// level of recursion.

        /// An alternate workaround is to compute both eigenvalues and
        /// eigenvectors, but that seems to be significantly slower.
        El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
        hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff = schur_complement_cholesky.blocks[block].Height() / 2 + 1;

        /// The default number of iterations is 40.  That is sometimes
        /// not enough, so we bump it up significantly.
        hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 400;

        El::DistMatrix<El::BigFloat, El::VR, El::STAR> eigenvalues(schur_complement_cholesky.blocks[block].Grid());
        El::HermitianEig(El::UpperOrLowerNS::LOWER, schur_complement_cholesky.blocks[block], eigenvalues,
                         hermitian_eig_ctrl);

        El::BigFloat max(0), min(std::numeric_limits<double>::max());
        for(int64_t row(0); row!=eigenvalues.Height(); ++row)
          {
            max=std::max(max,El::Abs(eigenvalues.Get(row,0)));
            min=std::min(min,El::Abs(eigenvalues.Get(row,0)));
          }
        El::BigFloat condition(max/min);
        if(condition>S_condition)
          {
            S_condition=condition;
            S_min=min;
            S_max=max;
          }
        // El::Print(schur_complement.blocks[block],"S");
        // std::cout << "\n";
        // exit(0);
        schur_complement_cholesky.blocks[block] = schur_complement.blocks[block];
      }
      

      Cholesky(El::UpperOrLowerNS::LOWER,
               schur_complement_cholesky.blocks[block]);
      cholesky_timer.stop();

      // SchurOffDiagonal = L^{-1} FreeVarMatrix
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
      //   Q = (L^{-1} B)^T (L^{-1} B)

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

  {
    std::stringstream ss;
    ss << "S Condition: "
       << El::mpi::Rank() << " "
       << S_condition << " "
       << S_max << " "
       << S_min << "\n";
    std::cout << ss.str();
  }
}
