#include "catch2/catch_amalgamated.hpp"

#include "test_util/test_util.hxx"
#include "unit_tests/util/util.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/Fmpz_Matrix.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/Matrix_Normalizer.hxx"

#include <vector>
#include <El.hpp>

using Test_Util::REQUIRE_Equal::diff;

std::vector<Blas_Job>
create_blas_jobs_split_remaining_primes(size_t num_ranks, size_t num_primes,
                                        El::Int output_matrix_height,
                                        size_t split_factor);

// Helper functions for calculating Q = P^T P
// using different methods
namespace
{
  // Calculate Q = P^T P using El::Syrk
  El::Matrix<El::BigFloat>
  calculate_matrix_square_El_Syrk(const El::Matrix<El::BigFloat> &P_matrix)
  {
    El::Matrix<El::BigFloat> Q_result_El_Syrk;
    El::Syrk(El::UpperOrLowerNS::UPPER, El::OrientationNS::TRANSPOSE,
             El::BigFloat(1.0), P_matrix, Q_result_El_Syrk);
    El::MakeSymmetric(El::UpperOrLowerNS::UPPER, Q_result_El_Syrk);
    return Q_result_El_Syrk;
  }

  // Calculate Q = P^T P using El::Gemm
  El::Matrix<El::BigFloat>
  calculate_matrix_square_El_Gemm(const El::Matrix<El::BigFloat> &P_matrix)
  {
    El::Matrix<El::BigFloat> Q_result_El_Gemm;
    El::Gemm(El::OrientationNS::TRANSPOSE, El::OrientationNS::NORMAL,
             El::BigFloat(1.0), P_matrix, P_matrix, Q_result_El_Gemm);
    return Q_result_El_Gemm;
  }

  // Calculate Q = P^T P using fmpz_mat_mul_blas:
  // - Normalize, convert to fmpz (big integers)
  // - Calculate residues (mod a bunch of primes)
  // - Multiply matrices of residues via BLAS
  // - Restore the result using Chinese Remainder Theorem
  // - Convert back to BigFloat, restore
  El::Matrix<El::BigFloat> calculate_matrix_square_fmpz_mat_mul_blas(
    const El::Matrix<El::BigFloat> &P_matrix, int bits)
  {
    auto P = P_matrix;
    INFO("auto P = P_matrix;");
    CAPTURE(P);
    CAPTURE(bits);
    // matrix only on this rank, thus El::mpi::COMM_SELF, no communication
    Matrix_Normalizer normalizer(P, bits, El::mpi::COMM_SELF);

    normalizer.normalize_and_shift_P(P);

    El::Matrix<El::BigFloat> PT;
    Transpose(P, PT);
    INFO("Transpose(P, PT);");

    Fmpz_Matrix PT_bigint(PT);
    Fmpz_Matrix P_bigint(P);
    Fmpz_Matrix Q(P.Width(), P.Width());
    fmpz_mat_mul_blas(Q.fmpz_matrix, PT_bigint.fmpz_matrix,
                      P_bigint.fmpz_matrix);

    El::Matrix<El::BigFloat> Q_result;
    Q.ToBigFloatMatrix(Q_result);

    normalizer.restore_Q(El::UPPER, Q_result);
    El::MakeSymmetric(El::UPPER, Q_result);
    return Q_result;
  }
}

TEST_CASE("calculate_Block_Matrix_square")
{
  INFO("input: dense tall NxK matrix P, splitted horizontally into blocks");
  INFO("output: NxN matrix Q := P^T P");
  INFO("We calculate Q with different methods, including our bigint_syrk_blas,"
       ", and compare the results.");

  MPI_Comm comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &comm);
  {
    INFO("TODO: the test currently doesn't work for several nodes");
    REQUIRE(El::mpi::Congruent(comm, El::mpi::COMM_WORLD));
  }

  CAPTURE(El::mpi::Rank());
  CAPTURE(El::mpi::Rank(comm));
  const El::Grid comm_grid(comm);

  int block_width = GENERATE(1, 10);
  DYNAMIC_SECTION("P_width=" << block_width)
  {
    int total_block_height = GENERATE(1, 10, 100, 1000);
    auto block_heights = Test_Util::random_split(total_block_height);
    size_t num_blocks = block_heights.size();

    DYNAMIC_SECTION("P_height=" << total_block_height
                                << " num_blocks=" << num_blocks)
    {
      bool use_dist_blocks = GENERATE(false, true);
      if(use_dist_blocks)
        INFO("Each block is a DistMatrix over a communicator");
      else
        INFO("Each block is stored on a single rank");

      DYNAMIC_SECTION("use_dist_blocks=" << use_dist_blocks)
      {
        int bits;
        CAPTURE(bits = El::gmp::Precision());
        int diff_precision;
        // In general, different calculations may give different results.
        // We cannot control the resulting precision exactly,
        // e.g. due to possible catastrophic cancellation.
        // If we are doing everything correctly, then we can expect
        // that results coincide up to bits/2 in most of the cases, this feels robust enough.
        // If we made a mistake, then usually the result is completely different,
        // so the test will fail.
        CAPTURE(diff_precision = bits / 2);

        // P_matrix is a tall matrix of all blocks.
        // We initialize it on rank=0 and then copy to all ranks.
        El::Matrix<El::BigFloat> P_matrix(total_block_height, block_width);
        El::Matrix<El::BigFloat> Q_result_El_Syrk(block_width, block_width);

        {
          INFO("Initialize P_matrix and calculate local Q with builtin "
               "Elemental and FLINT methods");
          if(El::mpi::Rank() == 0)
            {
              P_matrix
                = Test_Util::random_matrix(total_block_height, block_width);

              // Fill with 1.0 - for easier debug:
              // El::Fill(P_matrix, El::BigFloat(1.0));

              Q_result_El_Syrk = calculate_matrix_square_El_Syrk(P_matrix);

              // Double-check our result with Gemm
              auto Q_result_El_Gemm
                = calculate_matrix_square_El_Gemm(P_matrix);
              DIFF(Q_result_El_Syrk, Q_result_El_Gemm);

              // Check fmpz_mat_mul_blas, which is mathematically
              // the same as our bigint_syrk_blas method.
              auto Q_result_fmpz_mat_mul_blas
                = calculate_matrix_square_fmpz_mat_mul_blas(P_matrix, bits);
              DIFF_PREC(Q_result_El_Syrk, Q_result_fmpz_mat_mul_blas,
                        diff_precision);
            }

          //
          El::mpi::Broadcast(P_matrix.Buffer(), P_matrix.MemorySize(), 0,
                             El::mpi::COMM_WORLD);
          // Send copies of Q_result_El_Syrk to all ranks, to make comparison easy
          // TODO use DistMatrix instead?
          El::mpi::Broadcast(Q_result_El_Syrk.Buffer(),
                             Q_result_El_Syrk.MemorySize(), 0,
                             El::mpi::COMM_WORLD);
        }

        std::vector<El::DistMatrix<El::BigFloat>> P_matrix_blocks;
        std::vector<size_t> block_indices;
        {
          INFO("Initialize P_matrix_blocks for FLINT+BLAS "
               "multiplication");

          // Each block is either distributed over all comm,
          // or belongs to a single rank
          const auto &block_grid
            = use_dist_blocks ? comm_grid : El::Grid::Trivial();

          int global_block_offset = 0;
          for(size_t block_index = 0; block_index < block_heights.size();
              ++block_index)
            {
              int block_height = block_heights.at(block_index);

              El::DistMatrix<El::BigFloat> block(block_height, block_width,
                                                 block_grid);

              // If block goes for a single rank, then we assign blocks to ranks in a round-robin manner.
              if(use_dist_blocks
                 || block_index % El::mpi::Size(comm) == El::mpi::Rank(comm))
                {
                  for(int iLoc = 0; iLoc < block.LocalHeight(); ++iLoc)
                    for(int jLoc = 0; jLoc < block.LocalWidth(); ++jLoc)
                      {
                        int global_row
                          = block.GlobalRow(iLoc) + global_block_offset;
                        int global_col = block.GlobalCol(jLoc);
                        block.SetLocal(iLoc, jLoc,
                                       P_matrix.Get(global_row, global_col));
                      }
                  // TODO block indices for window!
                  block_indices.push_back(block_index);
                  P_matrix_blocks.push_back(block);
                }
              global_block_offset += block_height;
            }
        }

        // calculate Q = P^T P using bigint_syrk_blas
        // P is split into (split_factor) vertical bands
        for(size_t split_factor = 1; split_factor <= block_width;
            split_factor += 3)
          DYNAMIC_SECTION("split_factor=" << split_factor)
          {
            INFO("P matrix is split into " << split_factor
                                           << " vertical bands P_I");
            INFO("The blocks Q_IJ = P_I^T * P_J are calculated in parallel.");

            El::DistMatrix<El::BigFloat> Q_result(block_width, block_width,
                                                  comm_grid);
            {
              // blocks are distributed among all ranks, thus COMM_WORLD
              Matrix_Normalizer normalizer(P_matrix_blocks, block_width, bits,
                                           El::mpi::COMM_WORLD);
              CAPTURE(normalizer.column_norms);
              for(auto &block : P_matrix_blocks)
                {
                  normalizer.normalize_and_shift_P(block);
                }

              El::UpperOrLower uplo = El::UpperOrLowerNS::UPPER;

              auto create_jobs
                = [&split_factor](size_t num_ranks, size_t num_primes,
                                  int output_width) {
                    return create_blas_jobs_split_remaining_primes(
                      num_ranks, num_primes, output_width, split_factor);
                  };
              BigInt_Shared_Memory_Syrk_Context context(
                comm, bits, block_heights, block_width, block_indices,
                create_jobs);

              Timers timers(false);
              context.bigint_syrk_blas(uplo, P_matrix_blocks, Q_result,
                                       timers);
              {
                INFO("Check that normalized Q_ii = 1:");
                for(int iLoc = 0; iLoc < Q_result.LocalHeight(); ++iLoc)
                  for(int jLoc = 0; jLoc < Q_result.LocalWidth(); ++jLoc)
                    {
                      int i = Q_result.GlobalRow(iLoc);
                      int j = Q_result.GlobalCol(jLoc);
                      if(i == j)
                        {
                          CAPTURE(i);
                          auto value = Q_result.GetLocal(iLoc, jLoc);
                          DIFF_PREC(value >> 2 * normalizer.precision,
                                    El::BigFloat(1), diff_precision);
                        }
                    }
              }

              normalizer.restore_Q(uplo, Q_result);
              El::MakeSymmetric(uplo, Q_result);
            }

            // Check the result

            CAPTURE(Q_result);

            for(int iLoc = 0; iLoc < Q_result.LocalHeight(); ++iLoc)
              for(int jLoc = 0; jLoc < Q_result.LocalWidth(); ++jLoc)
                {
                  CAPTURE(iLoc);
                  CAPTURE(jLoc);
                  auto global_row = Q_result.GlobalRow(iLoc);
                  auto global_col = Q_result.GlobalCol(jLoc);
                  CAPTURE(global_row);
                  CAPTURE(global_col);

                  DIFF_PREC(Q_result.GetLocal(iLoc, jLoc),
                            Q_result_El_Syrk.Get(global_row, global_col),
                            diff_precision);
                }
          }
      }
    }
  }
}