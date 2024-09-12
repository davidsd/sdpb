#include "spectrum/Zeros.hxx"
#include "compute_lambda.hxx"
#include "pmp/PMP_Info.hxx"

#include <vector>

std::vector<El::BigFloat>
find_zeros(const El::Matrix<El::BigFloat> &c_minus_By_block,
           const PVM_Info &pvm, const El::BigFloat &threshold);

std::vector<Zeros>
compute_spectrum(const PMP_Info &pmp,
                 const std::vector<El::Matrix<El::BigFloat>> &c_minus_By,
                 const std::optional<std::vector<El::Matrix<El::BigFloat>>> &x,
                 const El::BigFloat &threshold, const bool &need_lambda)
{
  std::vector<Zeros> spectrum_blocks(pmp.blocks.size());
  for(size_t local_block_index = 0; local_block_index < pmp.blocks.size();
      ++local_block_index)
    {
      const auto &pvm_info = pmp.blocks.at(local_block_index);
      const auto &c_minus_By_block = c_minus_By.at(local_block_index);
      auto &spectrum_block = spectrum_blocks.at(local_block_index);

      spectrum_block.block_path = pvm_info.block_path;

      const auto zero_values
        = find_zeros(c_minus_By_block, pvm_info, threshold);
      if(need_lambda)
        {
          ASSERT(x.has_value());
          compute_lambda(
            pvm_info.sample_points, pvm_info.reduced_sample_scalings,
            pvm_info.reduced_prefactor, pvm_info.dim, x->at(local_block_index),
            zero_values, pvm_info.block_index, spectrum_block.zeros,
            spectrum_blocks.at(local_block_index).error);
        }
      else
        {
          for(auto &zero_value : zero_values)
            {
              spectrum_block.zeros.emplace_back(zero_value);
            }
        }
    }
  return spectrum_blocks;
}
