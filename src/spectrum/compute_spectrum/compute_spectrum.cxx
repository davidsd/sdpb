#include "spectrum/Zeros.hxx"
#include "compute_lambda.hxx"
#include "pmp/Polynomial_Matrix_Program.hxx"

#include <vector>

std::vector<El::BigFloat>
find_zeros(const El::Matrix<El::BigFloat> &c_minus_By_block,
           const Polynomial_Vector_Matrix &pvm, const El::BigFloat &threshold);

std::vector<Zeros>
compute_spectrum(const Polynomial_Matrix_Program &pmp,
                 const std::vector<El::Matrix<El::BigFloat>> &c_minus_By,
                 const std::optional<std::vector<El::Matrix<El::BigFloat>>> &x,
                 const El::BigFloat &threshold, const bool &need_lambda)
{
  std::vector<Zeros> spectrum_blocks(pmp.matrices.size());
  for(size_t block_index = 0; block_index < pmp.matrices.size(); ++block_index)
    {
      const auto &pvm = pmp.matrices.at(block_index);
      const auto &c_minus_By_block = c_minus_By.at(block_index);
      auto &spectrum_block = spectrum_blocks.at(block_index);

      spectrum_block.block_path = pmp.block_paths.at(block_index);

      const auto zero_values = find_zeros(c_minus_By_block, pvm, threshold);
      if(need_lambda)
        {
          ASSERT(x.has_value());
          compute_lambda(
            pvm.sample_points, pvm.reduced_sample_scalings,
            pvm.preconditioning_vector, pvm.reduced_prefactor,
            pvm.polynomials.Height(), x->at(block_index), zero_values,
            pmp.matrix_index_local_to_global.at(block_index),
            spectrum_block.zeros, spectrum_blocks.at(block_index).error);
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
