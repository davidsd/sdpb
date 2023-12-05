#include "sdp_read/sdp_read.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdp_convert/sdp_convert.hxx"

namespace fs = std::filesystem;

std::vector<El::Matrix<El::BigFloat>>
read_x(const fs::path &solution_dir,
       const std::vector<Polynomial_Vector_Matrix> &matrices)
{
  std::vector<El::Matrix<El::BigFloat>> result;
  result.reserve(matrices.size());
  const size_t num_procs(El::mpi::Size(El::mpi::COMM_WORLD));
  size_t x_index(El::mpi::Rank());
  for(auto &m : matrices)
    {
      result.emplace_back(m.sample_points.size() * m.rows * (m.rows + 1) / 2,
                          1);
      read_text_block(result.back(), solution_dir, "x_", x_index);
      x_index += num_procs;
    }
  return result;
}

std::vector<El::Matrix<El::BigFloat>>
read_x(const fs::path &solution_path,
       const std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  std::vector<El::Matrix<El::BigFloat>> result;
  result.reserve(matrices.size());
  const size_t num_procs(El::mpi::Size(El::mpi::COMM_WORLD));
  size_t x_index(El::mpi::Rank());
  for(auto &m : matrices)
    {
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
      read_text_block(result.back(), solution_path, "x_", x_index);
      x_index += num_procs;
    }
  return result;
}
