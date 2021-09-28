#include "../../../../SDP.hxx"
#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../Timers.hxx"
#include "../../../../../Boost_Float.hxx"
#include "../../../../../ostream_vector.hxx"

#include <apfp/config.h>
#include <apfp/matmul_library.h>

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

  static int64_t num_tries(0);
  const int64_t fpga_tries(1);
  for(size_t block = 0; block < schur_complement_cholesky.blocks.size();
      ++block)
    {
      auto &cholesky_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.cholesky_"
        + std::to_string(block_info.block_indices[block])));
      schur_complement_cholesky.blocks[block] = schur_complement.blocks[block];

      Cholesky(El::UpperOrLowerNS::LOWER,
               schur_complement_cholesky.blocks[block]);
      cholesky_timer.stop();

      // schur_off_diagonal = L^{-1} B
      auto &solve_timer(timers.add_and_start(
        "run.step.initializeSchurComplementSolver.Q.solve_"
        + std::to_string(block_info.block_indices[block])));

      schur_off_diagonal.blocks.push_back(sdp.free_var_matrix.blocks[block]);
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), schur_complement_cholesky.blocks[block],
               schur_off_diagonal.blocks[block]);

      solve_timer.stop();

      const size_t fpga_block(0);
      if(El::mpi::Rank()==0 && block==fpga_block && num_tries==fpga_tries)
      {
        std::string kernel_path(apfp::kBitstreamDirectory + std::string("/matmul_dram_kernel_hw.xclbin"));
        const int n_raw(schur_off_diagonal.blocks[block].Width()),
          k_raw(schur_off_diagonal.blocks[block].Height());

        const int n(n_raw + (apfp::kTileSize - n_raw % apfp::kTileSize) % apfp::kTileSize),
          k(k_raw + (apfp::kTileSize - k_raw % apfp::kTileSize) % apfp::kTileSize),
          m(n);

        mpfr_set_default_prec(apfp::kBits);
        mpfr_t zero;
        mpfr_init_set_ui(zero, 0, MPFR_RNDN);
        std::vector<mpfr_t> LB(k * n), LB_T(n*k), result(n*m);
        for(auto &element: LB)
          {
            mpfr_init_set(element, zero, MPFR_RNDN);
          }
        for(auto &element: LB_T)
          {
            mpfr_init_set(element, zero, MPFR_RNDN);
          }
        for(auto &element: result)
          {
            mpfr_init_set(element, zero, MPFR_RNDN);
          }
        std::stringstream ss;
        set_stream_precision(ss);
        El::BigFloat big_zero(0);
        size_t num_zeros(0);
        for(int64_t row(0); row != schur_off_diagonal.blocks[block].Height();
            ++row)
          for(int64_t column(0);
              column != schur_off_diagonal.blocks[block].Width(); ++column)
            {
              if(schur_off_diagonal.blocks[block].Get(row, column)==big_zero)
                {
                  ++num_zeros;
                }
              ss.str("");
              ss << schur_off_diagonal.blocks[block].Get(row, column);
              mpfr_init_set_str(LB_T[column + n * row], ss.str().c_str(), 10, MPFR_RNDN);
            }

        std::cout << "FPGA: " << block << " "
                  << n_raw << " "
                  << k_raw << " "
                  << n << " "
                  << k << " "
                  << m << " "
                  << n_raw*k_raw << " "
                  << num_zeros << " "
                  << (double(num_zeros)/(n_raw*k_raw)) << " "
                  << "\n" << std::flush;
        
        Timer fpga_total;
        apfp::Initialize(kernel_path, n, k, m);
        Timer fpga_copy_multiply;
        Timer fpga_copya;
        apfp::CopyAToDevice(LB_T.data());
        fpga_copya.stop();
        Timer fpga_copyb;
        apfp::CopyBToDevice(LB_T.data());
        fpga_copyb.stop();
        Timer fpga_multiply;
        apfp::MultiplyAxB();
        fpga_multiply.stop();
        Timer fpga_copyc;
        apfp::CopyCToHost(result.data());
        fpga_copyc.stop();
        fpga_copy_multiply.stop();
        fpga_total.stop();
        
        std::cout << "fpga: "
                  << fpga_total << " "
                  << fpga_copy_multiply << " "
                  << fpga_copya << " "
                  << fpga_copyb << " "
                  << fpga_multiply << " "
                  << fpga_copyc << " "
                  << "\n";
          
        // for(int64_t row(0); row != n_raw; ++row)
        //   for(int64_t column(0); column != n_raw; ++column)
        {
          // const int64_t row(12), column(17);
          // printf("result: %ld %ld ",row, column);
          // mpfr_out_str(stdout, 10, 0, result[row + n*column], MPFR_RNDD);
          // printf("\n");
        }

        apfp::Finalize();
      }
      // Q = (L^{-1} B)^T (L^{-1} B) = schur_off_diagonal^T schur_off_diagonal
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

      if(El::mpi::Rank()==0 && block==fpga_block && num_tries==fpga_tries)
        {
          std::cout << "Elemental: "
                    << syrk_timer << " "
                    // << Q_group_view.Get(12,17)
                    << "\n";
        }
    }
  ++num_tries;
}
