#include "diff.hxx"

#include "Parse_SDP.hxx"

#include <catch2/catch_amalgamated.hpp>

namespace fs = std::filesystem;

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void diff_sdp(const fs::path &a_sdp, const fs::path &b_sdp,
                unsigned int input_precision, unsigned int diff_precision,
                Test_Case_Runner runner, bool check_normalization)
  {
    INFO("diff sdp directories/zip files");
    REQUIRE(a_sdp != b_sdp);
    CAPTURE(a_sdp);
    CAPTURE(b_sdp);
    Float_Binary_Precision prec(input_precision, diff_precision);

    Parse_SDP a(a_sdp, runner);
    Parse_SDP b(b_sdp, runner);

    CAPTURE(a.sdp_dir);
    CAPTURE(b.sdp_dir);

    DIFF(a.control.num_blocks, b.control.num_blocks);
    // ignore "command", since it's unimportant:
    // diff(a.control.command, b.control.command);

    DIFF(a.objectives.constant, b.objectives.constant);
    DIFF(a.objectives.b, b.objectives.b);

    if(check_normalization)
      {
        DIFF(a.normalization.has_value(), b.normalization.has_value());
        if(a.normalization.has_value())
          DIFF(a.normalization->normalization, b.normalization->normalization);
      }

    DIFF(a.pmp_info.pmp_info, b.pmp_info.pmp_info);

    const auto num_blocks = a.control.num_blocks;
    {
      REQUIRE(a.block_info.size() == num_blocks);
      REQUIRE(b.block_info.size() == num_blocks);
      REQUIRE(a.block_data.size() == num_blocks);
      REQUIRE(b.block_data.size() == num_blocks);
    }

    for(size_t block_index = 0; block_index < num_blocks; ++block_index)
      {
        CAPTURE(block_index);

        const auto &a_block_info = a.block_info[block_index];
        const auto &b_block_info = b.block_info[block_index];
        DIFF(a_block_info.dim, b_block_info.dim);
        DIFF(a_block_info.num_points, b_block_info.num_points);

        const auto &a_block_data = a.block_data[block_index];
        const auto &b_block_data = b.block_data[block_index];
        DIFF(a_block_data.bilinear_bases_even,
             b_block_data.bilinear_bases_even);
        DIFF(a_block_data.bilinear_bases_odd, b_block_data.bilinear_bases_odd);
        DIFF(a_block_data.constraint_constants,
             b_block_data.constraint_constants);
        DIFF(a_block_data.constraint_matrix, b_block_data.constraint_matrix);
      }
  }
}
