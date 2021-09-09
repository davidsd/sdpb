#include "../sdp_read.hxx"
#include "../sdp_solve.hxx"
#include "../sdp_convert.hxx"

std::vector<El::Matrix<El::BigFloat>>
read_x(const boost::filesystem::path &solution_dir,
       const std::vector<Polynomial_Vector_Matrix> &matrices)
{
  std::vector<El::Matrix<El::BigFloat>> result;
  result.reserve(matrices.size());
  for(size_t block_index(0); block_index != matrices.size(); ++block_index)
    {
      auto &m(matrices[block_index]);
      result.emplace_back(m.sample_points.size() * m.rows * (m.rows + 1) / 2,
                          1);
      read_text_block(result.back(), solution_dir, "x_", block_index);
    }
  return result;
}

std::vector<El::Matrix<El::BigFloat>>
read_x(const boost::filesystem::path &solution_path,
       const std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  std::vector<El::Matrix<El::BigFloat>> result;
  result.reserve(matrices.size());
  for(size_t block_index(0); block_index != matrices.size(); ++block_index)
    {
      auto &m(matrices[block_index]);

      // Copied from sdp2input/write_output/write_output.cxx
      const size_t max_degree([&]() {
        int64_t result(0);
        for(auto &pvv : m.polynomials)
          for(auto &pv : pvv)
            for(auto &polynomial : pv)
              {
                result = std::max(result, polynomial.degree());
              }
        return result;
      }());

      result.emplace_back((max_degree + 1) * m.polynomials.size()
                            * (m.polynomials.size() + 1) / 2,
                          1);
      read_text_block(result.back(), solution_path, "x_", block_index);
    }
  return result;
}
