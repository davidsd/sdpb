#include "sdp_solve/sdp_solve.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_timing(const fs::path &checkpoint_out, const Block_Info &block_info,
                  const Timers &timers, const bool &debug,
                  El::Matrix<int32_t> &block_timings)
{
  if(debug)
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

  El::Zero(block_timings);
  std::string Q_prefix
    = "sdpb.solve.run.iter_2.step.initializeSchurComplementSolver.Q.";

  for(auto timer_prefix :
      {Q_prefix + "syrk_", Q_prefix + "solve_", Q_prefix + "cholesky_"})
    {
      for(auto &index : block_info.block_indices)
        {
          block_timings(index, 0) += timers.elapsed_milliseconds(
            timer_prefix + std::to_string(index));
        }
    }
  El::AllReduce(block_timings, El::mpi::COMM_WORLD);
  if(El::mpi::Rank() == 0)
    {
      if(debug)
        {
          El::Print(block_timings, "block_timings, ms:", ", ");
          El::Output();
        }

      fs::create_directories(checkpoint_out);
      fs::path block_timings_path(checkpoint_out / "block_timings");
      // Copy old block_timings to block_timings.0
      if(fs::exists(block_timings_path))
        {
          for(size_t index = 0; index < std::numeric_limits<size_t>::max();
              ++index)
            {
              fs::path backup_path = block_timings_path;
              backup_path.replace_extension("." + std::to_string(index));
              if(!fs::exists(backup_path))
                {
                  El::Output("Copy old ", block_timings_path, " to ",
                             backup_path);
                  fs::copy(block_timings_path, backup_path);
                  break;
                }
            }
        }

      std::ofstream block_timings_file(block_timings_path);
      for(int64_t row = 0; row < block_timings.Height(); ++row)
        {
          block_timings_file << block_timings(row, 0) << "\n";
        }
      ASSERT(block_timings_file.good(),
             "Error when writing to: ", block_timings_path);
    }
}
