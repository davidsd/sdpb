#include "catch2/catch_amalgamated.hpp"

#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/reduce_scatter_DistMatrix.hxx"
#include "sdpb_util/copy_matrix.hxx"
#include "test_util/test_util.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>
#include <vector>

using Test_Util::REQUIRE_Equal::diff;

void deallocate_unused_half(El::DistMatrix<El::BigFloat> &matrix,
                            El::UpperOrLower uplo);

void diff_uplo(const El::DistMatrix<El::BigFloat> &a,
               const El::DistMatrix<El::BigFloat> &b,
               std::optional<El::UpperOrLower> uplo)
{
  if(uplo.has_value())
    {
      diff(a, b, uplo.value());
    }
  else
    {
      DIFF(a, b);
    }
}

TEST_CASE("reduce_scatter_DistMatrix")
{
  INFO("Test reduce-scatter for DistMatrix.");
  INFO(
    "Input: DistMatrices for each \"node\" (node has node_comm_size ranks)");
  INFO("Input: global DistMatrix.");
  El::mpi::Comm comm_world = El::mpi::COMM_WORLD;

  int node_comm_size = GENERATE(1, 2, 3);

  DYNAMIC_SECTION("node_comm_size=" << node_comm_size)
  {
    CAPTURE(comm_world.Size());
    // Each node should have the same number of ranks
    REQUIRE(comm_world.Size() % node_comm_size == 0);
    int num_nodes = std::ceil((double)comm_world.Size() / node_comm_size);
    // Create quasi-nodes to emulate multi-node computations
    El::mpi::Comm node_comm;
    MPI_Comm_split(comm_world.comm, El::mpi::Rank() / node_comm_size, 0,
                   &node_comm.comm);
    int node_index = El::mpi::Rank() / node_comm_size;
    REQUIRE(node_comm.Size() == node_comm_size);

    CAPTURE(comm_world.Rank());
    CAPTURE(node_comm.Rank());
    const El::Grid node_grid(node_comm);

    for(std::optional<El::UpperOrLower> uplo :
        std::vector<std::optional<El::UpperOrLower>>{std::nullopt, El::UPPER,
                                                     El::LOWER})
      {
        int height = GENERATE(1, 5, 100);
        int width = uplo.has_value() ? height : height + 4;

        std::string uplo_string = "null";
        if(uplo.has_value())
          uplo_string = uplo.value() == El::UPPER ? "UPPER" : "LOWER";
        DYNAMIC_SECTION("uplo=" << uplo_string << " height=" << height
                                << " width=" << width)
        {
          El::DistMatrix<El::BigFloat> node_DistMatrix(height, width,
                                                       node_grid);
          El::DistMatrix<El::BigFloat> result_DistMatrix(height, width);
          El::DistMatrix<El::BigFloat> expected_result_DistMatrix(height,
                                                                  width);

          {
            INFO("Initialize matrices");
            std::vector<El::Matrix<El::BigFloat>> node_matrices(num_nodes);
            El::Matrix<El::BigFloat> node_Matrix;

            El::Matrix<El::BigFloat> expected_result_Matrix(height, width);
            for(size_t i = 0; i < num_nodes; ++i)
              {
                node_matrices.at(i) = Test_Util::random_matrix(height, width);
                El::Broadcast(node_matrices.at(i), El::mpi::COMM_WORLD, 0);

                if(i == node_index)
                  node_Matrix = node_matrices.at(i);

                expected_result_Matrix += node_matrices.at(i);
              }

            copy_matrix_from_root(expected_result_Matrix,
                                  expected_result_DistMatrix,
                                  El::mpi::COMM_WORLD);

            copy_matrix_from_root(node_Matrix, node_DistMatrix, node_comm);
          }

          {
            INFO("reduce-scatter");
            Timers timers;
            reduce_scatter(result_DistMatrix, node_DistMatrix, timers, uplo);

            // CAPTURE(result_DistMatrix);
            // CAPTURE(expected_result_DistMatrix);
            diff(result_DistMatrix, expected_result_DistMatrix, uplo);
          }
        }
      }
  }
}