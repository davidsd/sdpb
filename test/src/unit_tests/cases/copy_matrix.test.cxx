#include <catch2/catch_amalgamated.hpp>
#include "sdpb_util/copy_matrix.hxx"
#include "unit_tests/util/util.hxx"

using Test_Util::REQUIRE_Equal::diff;

TEST_CASE("copy_matrix")
{
  SECTION("copy_matrix_from_root")
  {
    int height = GENERATE(1, 10);
    int width = GENERATE(1, 21);
    DYNAMIC_SECTION("height=" << height << " width=" << width)
    {
      // original_matrix will be the same for all ranks
      El::InitializeRandom(true);
      El::Matrix<El::BigFloat> original_matrix(height, width);
      if(El::mpi::Rank() == 0)
        original_matrix = Test_Util::random_matrix(height, width);
      El::mpi::Broadcast(original_matrix.Buffer(),
                         original_matrix.MemorySize(), 0, El::mpi::COMM_WORLD);

      // input initialized only on root
      El::Matrix<El::BigFloat> input;
      if(El::mpi::Rank() == 0)
        {
          input = original_matrix;
        }
      else
        {
          // Set wrong size, this should lead to errors
          // if we try to access it on non-root rank.
          input.Resize(height + 1, width - 1);
        }

      El::DistMatrix<El::BigFloat> output(height, width);
      copy_matrix_from_root(input, output, output.DistComm());

      // compare output with original_matrix
      for(int iLoc = 0; iLoc < output.LocalHeight(); ++iLoc)
        for(int jLoc = 0; jLoc < output.LocalWidth(); ++jLoc)
          {
            int i = output.GlobalRow(iLoc);
            int j = output.GlobalCol(jLoc);
            CAPTURE(iLoc);
            CAPTURE(jLoc);
            CAPTURE(i);
            CAPTURE(j);
            DIFF(output.GetLocal(iLoc, jLoc), original_matrix.Get(i, j));
          }
    }
  }
}
