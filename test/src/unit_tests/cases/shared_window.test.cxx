#include "catch2/catch_amalgamated.hpp"

#include "sdpb_util/Shared_Window_Array.hxx"
#include "test_util/diff.hxx"

#include "sdp_solve/SDP_Solver/run/bigint_syrk/Block_Residue_Matrices_Window.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/Residue_Matrices_Window.hxx"

#include <El.hpp>

#include <vector>

using Test_Util::REQUIRE_Equal::diff;

TEST_CASE("MPI_Shared_Window")
{
  El::mpi::Comm comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &comm.comm);

  SECTION("Shared_Window_Array")
  {
    size_t size = 10;
    Shared_Window_Array<double> array(comm, size);

    for(size_t i = 0; i < size; ++i)
      {
        if(El::mpi::Rank(comm) == i % El::mpi::Size(comm))
          array[i] = i;
      }

    array.Fence();
    for(size_t i = 0; i < size; ++i)
      {
        DIFF(array[i], (double)i);
      }
  }

  SECTION("Residue_Matrices_Window")
  {
    size_t num_primes = 10;
    size_t height = 21;
    size_t width = 32;

    std::vector<El::Matrix<double>> matrices(num_primes);
    for(size_t p = 0; p < num_primes; ++p)
      {
        auto &matrix = matrices.at(p);
        matrix.Resize(height, width);
        for(size_t i = 0; i < height; ++i)
          for(size_t j = 0; j < width; ++j)
            matrix(i, j) = p + 11 * i + 23 * j;
      }

    // BLAS output will be stored here
    Residue_Matrices_Window<double> window(comm, num_primes, height, width);
    for(size_t p = 0; p < num_primes; ++p)
      for(size_t i = 0; i < height; ++i)
        for(size_t j = 0; j < width; ++j)
          {
            if(El::mpi::Rank(comm) == (i + p + j) % El::mpi::Size(comm))
              window.residues.at(p)(i, j) = matrices.at(p)(i, j);
          }

    window.Fence();

    for(size_t p = 0; p < num_primes; ++p)
      DIFF(window.residues.at(p), matrices.at(p));
  }

  SECTION("Block_Residue_Matrices_Window")
  {
    size_t num_primes = 10;
    size_t num_blocks = 21;
    size_t width = 32;

    std::vector<std::vector<El::Matrix<double>>> block_residues(num_primes);
    std::vector<El::Int> block_heights(num_blocks);
    for(size_t b = 0; b < num_blocks; ++b)
      {
        block_heights.at(b) = width + b; // to make heights different
      }

    Block_Residue_Matrices_Window<double> window(comm, num_primes, num_blocks,
                                                 block_heights, width);
    for(size_t p = 0; p < num_primes; ++p)
      {
        block_residues.at(p).resize(num_blocks);
        for(size_t b = 0; b < num_blocks; ++b)
          {
            auto &matrix = block_residues.at(p).at(b);
            size_t height = block_heights.at(b);
            matrix.Resize(height, width);
            for(size_t i = 0; i < height; ++i)
              {
                for(size_t j = 0; j < width; ++j)
                  {
                    matrix(i, j) = -1.0 + p + b * 10 + i * 100 + j * 1001;
                    if(El::mpi::Rank(comm)
                       == (p + b + i + j) % El::mpi::Size(comm))
                      window.block_residues.at(p).at(b)(i, j) = matrix(i, j);
                  }
              }
          }
      }

    window.Fence();

    for(size_t p = 0; p < num_primes; ++p)
      {
        CAPTURE(p);
        size_t global_start_row = 0;
        for(size_t b = 0; b < num_blocks; ++b)
          {
            CAPTURE(b);
            DIFF(window.block_residues[p][b], block_residues[p][b]);

            INFO("Check that window.block_residues contains "
                 "submatrices of window.residues:");
            auto height = block_residues[p][b].Height();
            auto width = block_residues[p][b].Width();
            for(int i = 0; i < height; ++i)
              {
                CAPTURE(i);
                for(int j = 0; j < width; ++j)
                  {
                    CAPTURE(j);
                    DIFF(window.block_residues[p][b](i, j),
                         window.residues[p](i + global_start_row, j));
                  }
              }
            global_start_row += height;
          }
      }
  }
}
