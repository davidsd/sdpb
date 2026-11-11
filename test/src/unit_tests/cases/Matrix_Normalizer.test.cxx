#include "catch2/catch_amalgamated.hpp"

#include "sdp_solve/SDP_Solver/run/bigint_syrk/Matrix_Normalizer.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>

using Test_Util::REQUIRE_Equal::diff;

TEST_CASE("normalize_and_shift")
{
  int height = GENERATE(1, 10, 100);

  DYNAMIC_SECTION("height=" << height)
  {
    int bits = El::gmp::Precision();
    int diff_precision;
    CAPTURE(diff_precision = bits / 2);

    // width is not really important,
    // we just create non-square matrix for generality
    int width = height / 2 + 1;

    // non-square matrix for
    El::DistMatrix<El::BigFloat> P_matrix
      = Test_Util::random_distmatrix(height, width);
    CAPTURE(P_matrix.Height());
    CAPTURE(P_matrix.Width());
    CAPTURE(P_matrix.LocalHeight());
    CAPTURE(P_matrix.LocalWidth());

    Block_Matrix P_block_matrix({}, P_matrix.Width(), {}, P_matrix.Grid());

    // Split P_matrix horizontally
    // into blocks with different random heights
    auto block_heights = Test_Util::random_split(P_matrix.Height());
    size_t num_blocks = block_heights.size();
    CAPTURE(num_blocks);
    CAPTURE(block_heights);

    {
      INFO("Split P_matrix into horizontal bands");
      int begin_row = 0;
      for(size_t block_index = 0; block_index < num_blocks; ++block_index)
        {
          int block_height = block_heights.at(block_index);
          int end_row = begin_row + block_height;

          ASSERT(block_height > 0);
          ASSERT(end_row > begin_row);
          ASSERT(end_row <= P_matrix.Height());

          El::Range<int> rows(begin_row, end_row);
          El::Range<int> cols(0, P_matrix.Width());
          auto block = P_matrix(rows, cols);
          P_block_matrix.blocks.push_back(block);

          begin_row = end_row;
        }
      CAPTURE(P_block_matrix.blocks.size());
    }

    auto initial_P_matrix = P_matrix;

    const auto P_normalizer
      = normalize_and_shift<Matrix_Normalization_Kind::COLUMNS>(P_matrix,
                                                                bits);
    const auto P_block_normalizer
      = normalize_and_shift<Matrix_Normalization_Kind::COLUMNS>(
        P_block_matrix, bits, P_matrix.DistComm());
    CAPTURE(P_normalizer.bits);
    CAPTURE(P_normalizer.norms);
    CAPTURE(P_block_normalizer.norms);

    REQUIRE(P_matrix.Get(0, 0) != initial_P_matrix.Get(0, 0));
    DIFF(P_block_matrix.blocks.at(0).Get(0, 0), P_matrix.Get(0, 0));

    SECTION("restore_Q")
    {
      INFO(
        "Check that calculating Q with and without intermediate normalization "
        "gives the same result up to a reasonable diff_precision");

      El::DistMatrix<El::BigFloat> initial_Q, Q;
      El::UpperOrLower uplo = El::UPPER;
      El::Syrk(uplo, El::OrientationNS::TRANSPOSE, El::BigFloat(1),
               initial_P_matrix, initial_Q);
      El::Syrk(uplo, El::OrientationNS::TRANSPOSE, El::BigFloat(1), P_matrix,
               Q);

      {
        INFO("Check that normalized matrix squared has 1.0 on diagonal");
        for(int i = 0; i < P_matrix.LocalHeight(); ++i)
          for(int j = 0; j < P_matrix.LocalWidth(); ++j)
            {
              if(Q.GlobalRow(i) == Q.GlobalCol(j))
                DIFF_PREC(Q.GetLocal(i, j) >> 2 * bits, El::BigFloat(1.0),
                          diff_precision);
            }
      }

      auto Q_from_P = Q;
      auto Q_from_P_blocks = Q;

      El::MakeSymmetric(uplo, initial_Q);

      restore_syrk_output(uplo, P_normalizer, Q_from_P);
      El::MakeSymmetric(uplo, Q_from_P);
      DIFF_PREC(initial_Q, Q_from_P, diff_precision);

      restore_syrk_output(uplo, P_block_normalizer, Q_from_P_blocks);
      El::MakeSymmetric(uplo, Q_from_P_blocks);
      DIFF_PREC(initial_Q, Q_from_P_blocks, diff_precision);
    }

    SECTION("restore_P")
    {
      {
        INFO("restore_P");
        P_normalizer.restore(P_matrix);
        DIFF(P_matrix, initial_P_matrix);
      }

      {
        INFO("restore_P_blocks");
        P_block_normalizer.restore(P_block_matrix);

        int begin_row = 0;
        for(size_t block_index = 0; block_index < num_blocks; ++block_index)
          {
            int block_height = block_heights.at(block_index);
            int end_row = begin_row + block_height;

            El::Range<int> rows(begin_row, end_row);
            El::Range<int> cols(0, P_matrix.Width());
            auto block = P_matrix(rows, cols);

            CAPTURE(block_index);
            CAPTURE(begin_row);
            CAPTURE(block_height);
            DIFF(block, P_block_matrix.blocks.at(block_index));

            begin_row = end_row;
          }
      }
    }
  }
}