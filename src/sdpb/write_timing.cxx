#include "Block_Info.hxx"
#include "../Timers.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void write_timing(const boost::filesystem::path &checkpoint_out,
                  const Block_Info &block_info, const Timers &timers,
                  const bool &debug)
{
  if(debug)
    {
      timers.write_profile(checkpoint_out.string() + ".profiling."
                           + std::to_string(El::mpi::Rank()));
    }

  El::Matrix<int32_t> block_timings(block_info.dimensions.size(), 1);
  El::Zero(block_timings);
  for(auto &index : block_info.block_indices)
    {
      block_timings(index, 0) = timers.elapsed_milliseconds(
                                  "run.step.initializeSchurComplementSolver."
                                  "Qcomputation.syrk_"
                                  + std::to_string(index))
                                + timers.elapsed_milliseconds(
                                    "run.step.initializeSchurComplementSolver."
                                    "Qcomputation.solve_"
                                    + std::to_string(index))
                                + timers.elapsed_milliseconds(
                                    "run.step.initializeSchurComplementSolver."
                                    "Qcomputation.cholesky_"
                                    + std::to_string(index));
    }
  El::AllReduce(block_timings, El::mpi::COMM_WORLD);
  if(El::mpi::Rank() == 0)
    {
      boost::filesystem::create_directories(checkpoint_out);
      boost::filesystem::ofstream block_timings_file(checkpoint_out
                                                     / "block_timings");
      for(int64_t row = 0; row < block_timings.Height(); ++row)
        {
          block_timings_file << block_timings(row, 0) << "\n";
        }
    }
}
