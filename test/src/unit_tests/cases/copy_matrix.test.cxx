#include <catch2/catch_amalgamated.hpp>
#include "sdpb_util/copy_matrix.hxx"
#include "unit_tests/util/util.hxx"

using Test_Util::REQUIRE_Equal::diff;

void copy_matrix_from_root_impl_send_recv(
  const El::Matrix<El::BigFloat> &source,
  El::DistMatrix<El::BigFloat> &destination, const El::mpi::Comm &comm);

void copy_matrix_from_root_impl_shared_window(
  const El::Matrix<El::BigFloat> &source,
  El::DistMatrix<El::BigFloat> &destination, const El::mpi::Comm &comm);

TEST_CASE("copy_matrix")
{
  SECTION("copy_matrix_from_root")
  {
    int height = GENERATE(1, 10);
    int width = GENERATE(1, 21);
    bool use_shared_window = GENERATE(false, true);

    DYNAMIC_SECTION((use_shared_window ? "copy_matrix_from_root_shared_window"
                                       : "copy_matrix_from_root_send_recv"))
    DYNAMIC_SECTION("height=" << height << " width=" << width)
    {
      const auto &copy_matrix_func
        = use_shared_window ? copy_matrix_from_root_impl_shared_window
                            : copy_matrix_from_root_impl_send_recv;

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
        input = original_matrix;

      El::DistMatrix<El::BigFloat> output(height, width);
      copy_matrix_func(input, output, output.DistComm());

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
