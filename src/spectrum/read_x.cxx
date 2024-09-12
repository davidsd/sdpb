#include "pmp/PMP_Info.hxx"
#include "sdp_solve/read_text_block.hxx"
#include "sdpb_util/Timers/Timers.hxx"

namespace fs = std::filesystem;

std::vector<El::Matrix<El::BigFloat>>
read_x(const fs::path &solution_path, const PMP_Info &pmp, Timers &timers)
{
  Scoped_Timer timer(timers, "read_x");
  std::vector<El::Matrix<El::BigFloat>> result;
  result.reserve(pmp.blocks.size());
  for(auto &block : pmp.blocks)
    {
      const auto height = block.dim;
      const auto width = height;
      result.emplace_back(
        block.sample_points.size() * width * (height + 1) / 2, 1);
      read_text_block(result.back(), solution_path, "x_", block.block_index);
    }
  return result;
}
