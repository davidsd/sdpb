#include "sdp_solve/sdp_solve.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_profiling(const fs::path &checkpoint_out, const Timers &timers)
{
  // Write profiling for each rank to ck.profiling/profiling.{rank}
  fs::path parent_dir = checkpoint_out.string() + ".profiling";

  // Move old profiling data to ck.profiling.0/
  if(El::mpi::Rank() == 0 && fs::exists(parent_dir))
    {
      for(size_t index = 0; index < std::numeric_limits<size_t>::max();
          ++index)
        {
          const fs::path backup_dir
            = parent_dir.string() + "." + std::to_string(index);
          if(!fs::exists(backup_dir))
            {
              El::Output("Move old ", parent_dir, " to ", backup_dir);
              fs::rename(parent_dir, backup_dir);
              break;
            }
        }
    }
  // Barrier to ensure that we'll move old ck.profiling/ before writing to it
  El::mpi::Barrier();
  timers.write_profile(parent_dir
                       / ("profiling." + std::to_string(El::mpi::Rank())));
}

void write_block_timings(const fs::path &checkpoint_out,
                         const Block_Info &block_info,
                         const El::Matrix<int32_t> &block_timings_ms,
                         const Verbosity verbosity)
{
  if(El::mpi::Rank() != 0)
    return;

  if(verbosity >= Verbosity::debug)
    {
      El::Print(block_timings_ms, "block_timings, ms:", ", ");
      El::Output();
    }
  ASSERT_EQUAL(block_timings_ms.Height(), block_info.dimensions.size());
  ASSERT_EQUAL(block_timings_ms.Width(), 1);

  fs::create_directories(checkpoint_out);
  fs::path block_timings_path(checkpoint_out / "block_timings");
  if(fs::exists(block_timings_path))
    {
      PRINT_WARNING("block_timings file exists and will be overwritten: ",
                    block_timings_path);
    }

  std::ofstream block_timings_file(block_timings_path);
  for(int64_t row = 0; row < block_timings_ms.Height(); ++row)
    {
      const auto block_time = block_timings_ms(row, 0);
      ASSERT(block_time >= 0, "block = ", row, ": block_time = ", block_time,
             "ms is negative, probably overflow occurred!");
      block_timings_file << block_time << "\n";
    }
  ASSERT(block_timings_file.good(),
         "Error when writing to: ", block_timings_path);
}
