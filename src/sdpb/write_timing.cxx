#include "sdp_solve/sdp_solve.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_timing(const fs::path &checkpoint_out, const Block_Info &block_info,
                  const Timers &timers, const bool &debug,
                  const El::Matrix<int32_t> &block_timings_ms)
{
  if(debug)
    {
      timers.write_profile(checkpoint_out.string() + ".profiling/profiling."
                           + std::to_string(El::mpi::Rank()));
    }
  if(El::mpi::Rank() == 0)
    {
      if(debug)
        {
          El::Print(block_timings_ms, "block_timings, ms:", ", ");
          El::Output();
        }

      ASSERT_EQUAL(block_timings_ms.Height(), block_info.dimensions.size());
      ASSERT_EQUAL(block_timings_ms.Width(), 1);

      fs::create_directories(checkpoint_out);
      fs::path block_timings_path(checkpoint_out / "block_timings");
      std::ofstream block_timings_file(block_timings_path);
      for(int64_t row = 0; row < block_timings_ms.Height(); ++row)
        {
          const auto block_time = block_timings_ms(row, 0);
          ASSERT(block_time >= 0, "block = ", row,
                 ": block_time = ", block_time,
                 "ms is negative, probably overflow occurred!");
          block_timings_file << block_time << "\n";
        }
      ASSERT(block_timings_file.good(),
             "Error when writing to: ", block_timings_path);
    }
}
