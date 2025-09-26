#include "spectrum/Zeros.hxx"
#include "compute_lambda.hxx"
#include "pmp/PMP_Info.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>
#include <vector>

std::vector<El::BigFloat>
find_zeros(const El::Matrix<El::BigFloat> &c_minus_By_block,
           const PVM_Info &pvm, const El::BigFloat &threshold,
           const El::BigFloat &max_zero, Timers &timers);

void write_profiling(const std::filesystem::path &spectrum_output_path,
                     Timers &timers);

std::vector<Zeros>
compute_spectrum(const PMP_Info &pmp,
                 const std::vector<El::Matrix<El::BigFloat>> &c_minus_By,
                 const std::optional<std::vector<El::Matrix<El::BigFloat>>> &x,
                 const El::BigFloat &threshold, const El::BigFloat &max_zero,
                 const bool &need_lambda, const Verbosity &verbosity,
                 const std::filesystem::path &output_path, Timers &timers)
{
  Scoped_Timer timer(timers, "compute_spectrum");
  std::vector<Zeros> spectrum_blocks(pmp.blocks.size());
  for(size_t local_block_index = 0; local_block_index < pmp.blocks.size();
      ++local_block_index)
    {
      const auto &pvm_info = pmp.blocks.at(local_block_index);
      const auto block_index = pvm_info.block_index;

      Scoped_Timer block_timer(timers, "block_" + std::to_string(block_index));

      const auto &c_minus_By_block = c_minus_By.at(local_block_index);
      auto &spectrum_block = spectrum_blocks.at(local_block_index);

      spectrum_block.block_path = pvm_info.block_path;

      try
        {
          const auto zero_values = find_zeros(c_minus_By_block, pvm_info,
                                              threshold, max_zero, timers);
          if(need_lambda)
            {
              ASSERT(x.has_value());
              compute_lambda(pvm_info, x->at(local_block_index), zero_values,
                             spectrum_block, timers);
            }
          else
            {
              for(auto &zero_value : zero_values)
                {
                  spectrum_block.zeros.emplace_back(zero_value);
                }
            }
        }
      catch(const std::exception &e)
        {
          PRINT_WARNING("Failed to compute spectrum for block_",
                        pvm_info.block_index,
                        " block_path=", pvm_info.block_path, ":\n", e.what());
        }
      if(verbosity >= Verbosity::trace)
        {
          double time = block_timer.elapsed_milliseconds() / 1000.0;
          El::Output("Finished block_", block_index, ": dim=", pvm_info.dim,
                     " num_points=", pvm_info.sample_points.size(),
                     " num_zeros=", spectrum_block.zeros.size(),
                     " time=", time, "s");
          // Update profiling data.
          // Can be helpful to look at profiling for unfinished runs,
          // if some blocks take much longer than expected.
          write_profiling(output_path, timers);
        }
    }
  return spectrum_blocks;
}
