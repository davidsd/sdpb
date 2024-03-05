#include "pmp_read/pmp_read.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "pmp2sdp/write_sdp.hxx"

namespace fs = std::filesystem;

std::vector<El::Matrix<El::BigFloat>>
read_x(const fs::path &solution_path,
       const std::vector<Polynomial_Vector_Matrix> &matrices,
       const std::vector<size_t> &block_indices)
{
  std::vector<El::Matrix<El::BigFloat>> result;
  result.reserve(matrices.size());
  for(size_t i = 0; i < matrices.size(); ++i)
    {
      const auto block_index = block_indices.at(i);
      const auto &m = matrices.at(i);
      result.emplace_back(m.sample_points.size() * m.polynomials.Width() * (m.polynomials.Height() + 1) / 2,
                          1);
      read_text_block(result.back(), solution_path, "x_", block_index);
    }
  return result;
}
