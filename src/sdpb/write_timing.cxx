#include "Timers.hxx"
#include "Block_Info.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void write_timing(const boost::filesystem::path &out_file,
                  const boost::filesystem::path &block_timings_filename,
                  const Block_Info &block_info, const Timers &timers)
{
  timers.write_profile(out_file.string() + ".profiling."
                       + std::to_string(El::mpi::Rank()));
  if(!block_timings_filename.empty())
    {
      El::Matrix<int32_t> block_timings(block_info.dimensions.size(), 1);
      El::Zero(block_timings);
      for(auto &index : block_info.block_indices)
        {
          block_timings(index, 0)
            = timers.elapsed_milliseconds(
                "run.step.initializeSchurComplementSolver."
                "Qcomputation.syrk_"
                + std::to_string(index))
              + timers.elapsed_milliseconds(
                  "run.step.initializeSchurComplementSolver."
                  "Qcomputation.solve_"
                  + std::to_string(index));
        }
      El::AllReduce(block_timings, El::mpi::COMM_WORLD);
      if(El::mpi::Rank() == 0)
        {
          boost::filesystem::ofstream block_timings_file(
            block_timings_filename);
          for(int64_t row = 0; row < block_timings.Height(); ++row)
            {
              block_timings_file << block_timings(row, 0) << "\n";
            }
        }
    }
}
