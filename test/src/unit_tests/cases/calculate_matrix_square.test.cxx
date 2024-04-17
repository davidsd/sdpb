#include "catch2/catch_amalgamated.hpp"

#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/Matrix_Normalizer.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/fmpz/Fmpz_Matrix.hxx"
#include "test_util/test_util.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>
#include <vector>

using Test_Util::REQUIRE_Equal::diff;

Blas_Job_Schedule create_blas_job_schedule_split_remaining_primes(
  Blas_Job::Kind kind, El::UpperOrLower uplo, size_t num_ranks,
  size_t num_primes, El::Int output_matrix_height, El::Int output_matrix_width,
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
    if(fmpz_mat_mul_blas(Q.fmpz_matrix, PT_bigint.fmpz_matrix,
                         P_bigint.fmpz_matrix)
       == 0)
      {
        INFO(
          "fmpz_mat_mul_blas() returned 0 (probably since FLINT was compiled "
          "without BLAS), falling back to fmpz_mat_mul_multi_mod()");
        fmpz_mat_mul_multi_mod(Q.fmpz_matrix, PT_bigint.fmpz_matrix,
                               P_bigint.fmpz_matrix);
      }

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
       " and compare the results.");

  El::mpi::Comm comm_world = El::mpi::COMM_WORLD;
  CAPTURE(comm_world.Size());

  {
    El::mpi::Comm comm_shmem;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &comm_shmem.comm);
    {
      INFO("TODO: the test currently doesn't work for several nodes");
      REQUIRE(El::mpi::Congruent(comm_shmem, comm_world));
    }
  }

  INFO("NB: you should run the test on 6/12/18/... ranks,"
       " e.g. mpirun -n 6 unit_tests");
  const int node_size = GENERATE(1, 2, 3, 6);
  for(int group_size : std::set{1, node_size})
    {
      DYNAMIC_SECTION("node_size=" << node_size
                                   << " group_size=" << group_size)
      {
        {
          INFO("Total number of ranks should be multiple of node_size.");
          REQUIRE(comm_world.Size() % node_size == 0);
        }
        int num_nodes = comm_world.Size() / node_size;
        CAPTURE(num_nodes);
        // Create quasi-nodes to emulate multi-node computations
        El::mpi::Comm node_comm;
        MPI_Comm_split(comm_world.comm, El::mpi::Rank() / node_size, 0,
                       &node_comm.comm);
        int node_index = comm_world.Rank() / node_size;
        REQUIRE(node_comm.Size() == node_size);

        CAPTURE(comm_world.Rank());
        CAPTURE(node_comm.Rank());
        const El::Grid node_grid(node_comm);

        El::mpi::Comm group_comm;
        MPI_Comm_split(node_comm.comm, node_comm.Rank() / group_size, 0,
                       &group_comm.comm);
        CAPTURE(group_comm.Rank());
        const int group_index_in_node = node_comm.Rank() / group_size;
        REQUIRE(node_comm.Size() % group_size == 0);
        const int num_groups_per_node = node_comm.Size() / group_size;
        const int num_groups_global = num_groups_per_node * num_nodes;
        const El::Grid group_grid(group_comm);

        int block_width = GENERATE(1, 10);
        DYNAMIC_SECTION("P_width=" << block_width)
        {
          int total_block_height = GENERATE(1, 10, 100, 1000);
          auto block_heights = Test_Util::random_split(total_block_height);
          if(block_heights.size() < num_groups_global)
            {
              // Each group should have at least one block
              total_block_height += num_groups_global - block_heights.size();
              block_heights.resize(num_groups_global, 1);
            }
          CAPTURE(block_heights);
          size_t num_blocks = block_heights.size();

          auto num_primes
            = Fmpz_Comb(El::gmp::Precision(), El::gmp::Precision(), 1,
                        total_block_height)
                .num_primes;

          size_t max_shared_memory_bytes
            = GENERATE(std::numeric_limits<size_t>::max(), 1, 0);
          // We cannot use variables inside GENERATE, e.g. max_shared_memory_bytes = GENERATE(block_width)
          // Thus we set them below:
          INFO((max_shared_memory_bytes == 0
                  ? "Splitting both P and Q memory windows"
                : max_shared_memory_bytes == 1
                  ? "Splitting P memory window"
                  : "Do not split memory windows"));
          if(max_shared_memory_bytes == 1)
            {
              // Do not split output window
              size_t output_window_height = block_width;
              size_t output_window_width = output_window_height;
              // Input window will have 3 rows per MPI group
              size_t input_window_height = 3 * El::mpi::Size();
              size_t input_window_width = block_width;
              max_shared_memory_bytes
                = (output_window_height * output_window_width
                   + input_window_height * input_window_width)
                  * num_primes * sizeof(double);
            }
          else if(max_shared_memory_bytes == 0)
            {
              // Split output window by 3x3 submatrices
              size_t output_window_height = std::ceil(block_width / 3.0);
              size_t output_window_width = output_window_height;
              // Two identical input windows of minimal size, i.e. one row per MPI group
              size_t input_window_height = El::mpi::Size();
              size_t input_window_width = output_window_width;
              max_shared_memory_bytes
                = (output_window_height * output_window_width
                   + 2 * input_window_height * input_window_width)
                  * num_primes * sizeof(double);
            }

          DYNAMIC_SECTION("P_height="
                          << total_block_height << " num_blocks=" << num_blocks
                          << " maxSharedMemory=" << max_shared_memory_bytes)
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
            El::Matrix<El::BigFloat> Q_result_El_Syrk(block_width,
                                                      block_width);

            {
              INFO("Initialize P_matrix and calculate local Q with builtin "
                   "Elemental and FLINT methods");
              if(El::mpi::Rank() == 0)
                {
                  P_matrix = Test_Util::random_matrix(total_block_height,
                                                      block_width);

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
                    = calculate_matrix_square_fmpz_mat_mul_blas(P_matrix,
                                                                bits);
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

              int global_block_offset = 0;
              for(size_t block_index = 0; block_index < block_heights.size();
                  ++block_index)
                {
                  int block_height = block_heights.at(block_index);

                  El::DistMatrix<El::BigFloat> block(block_height, block_width,
                                                     group_grid);

                  // round-robin
                  if(block_index % num_nodes == node_index
                     && block_index % num_groups_per_node
                          == group_index_in_node)
                    {
                      for(int iLoc = 0; iLoc < block.LocalHeight(); ++iLoc)
                        for(int jLoc = 0; jLoc < block.LocalWidth(); ++jLoc)
                          {
                            int global_row
                              = block.GlobalRow(iLoc) + global_block_offset;
                            int global_col = block.GlobalCol(jLoc);
                            block.SetLocal(
                              iLoc, jLoc,
                              P_matrix.Get(global_row, global_col));
                          }
                      // TODO block indices for window!
                      block_indices.push_back(block_index);
                      P_matrix_blocks.push_back(block);
                    }
                  global_block_offset += block_height;
                }
              // Each MPI group should get at least 1 block
              REQUIRE(!P_matrix_blocks.empty());
            }

            // calculate Q = P^T P using bigint_syrk_blas
            // P is split into (blas_schedule_split_factor) vertical bands
            for(size_t blas_schedule_split_factor = 1;
                blas_schedule_split_factor <= block_width;
                blas_schedule_split_factor += 3)
              DYNAMIC_SECTION(
                "blas_split_factor=" << blas_schedule_split_factor)
              {
                INFO("P matrix is split into " << blas_schedule_split_factor
                                               << " vertical bands P_I");
                INFO("The blocks Q_IJ = P_I^T * P_J are calculated in "
                     "parallel.");

                El::DistMatrix<El::BigFloat> Q_result(block_width,
                                                      block_width);
                {
                  // blocks are distributed among all ranks, thus COMM_WORLD
                  Matrix_Normalizer normalizer(P_matrix_blocks, block_width,
                                               bits, El::mpi::COMM_WORLD);
                  CAPTURE(normalizer.column_norms);
                  for(auto &block : P_matrix_blocks)
                    {
                      normalizer.normalize_and_shift_P(block);
                    }

                  const El::UpperOrLower uplo = El::UpperOrLowerNS::UPPER;

                  auto create_job_schedule
                    = [&blas_schedule_split_factor](
                        Blas_Job::Kind kind, El::UpperOrLower uplo,
                        size_t num_ranks, size_t num_primes, int output_height,
                        int output_width, Verbosity verbosity) {
                        return create_blas_job_schedule_split_remaining_primes(
                          kind, uplo, num_ranks, num_primes, output_height,
                          output_width, blas_schedule_split_factor);
                      };

                  // NB: group_comm.Size() is the same for everyone,
                  // otherwise we'd have to synchronize it
                  std::vector<int> group_comm_sizes_per_node(
                    num_groups_per_node, group_size);
                  std::vector<El::Int> blocks_height_per_group(
                    num_groups_per_node);
                  for(const auto &block : P_matrix_blocks)
                    {
                      if(block.DistComm().Rank() == 0)
                        blocks_height_per_group.at(group_index_in_node)
                          += block.Height();
                    }
                  El::mpi::AllReduce(blocks_height_per_group.data(),
                                     num_groups_per_node, node_comm);

                  {
                    const Verbosity verbosity = Verbosity::regular;
                    BigInt_Shared_Memory_Syrk_Context context(
                      node_comm, group_index_in_node,
                      group_comm_sizes_per_node, bits, max_shared_memory_bytes,
                      blocks_height_per_group, block_width, block_indices,
                      verbosity, create_job_schedule);

                    Timers timers;
                    El::Matrix<int32_t> block_timings_ms(num_blocks, 1);
                    El::Zero(block_timings_ms);
                    INFO("Calling bigint_syrk_blas()...");
                    context.bigint_syrk_blas(uplo, P_matrix_blocks, Q_result,
                                             timers, block_timings_ms);
                  }
                  {
                    INFO("Check that normalized Q_ii = 1:");
                    // TODO: if check fails only on some rank!=0,
                    // rank=0 will hang, waiting for
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
}