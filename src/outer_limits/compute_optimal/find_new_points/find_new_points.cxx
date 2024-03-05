#include "outer_limits/Function.hxx"
#include "sdpb_util/Mesh.hxx"

#include <El.hpp>

std::vector<El::BigFloat>
get_new_points(const Mesh &mesh, const El::BigFloat &block_epsilon);

El::BigFloat
eval_summed(const El::BigFloat &epsilon, const El::BigFloat &infinity,
            const std::vector<std::vector<Function>> &summed_functions,
            const El::BigFloat &x);

void find_new_points(
  const size_t &num_blocks, const size_t &rank, const size_t &num_procs,
  const El::BigFloat &mesh_threshold, const El::BigFloat &epsilon,
  const El::BigFloat &infinity,
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<El::BigFloat> &weights,
  const std::vector<std::set<El::BigFloat>> &points,
  std::vector<std::vector<El::BigFloat>> &new_points, bool &has_new_points)
{
  std::vector<size_t> num_new_points(num_blocks, 0);
  for(size_t block(rank); block < num_blocks; block += num_procs)
    {
      // TODO: These can both be precomputed
      El::BigFloat max_delta(infinity), block_scale(0);
      size_t max_degree(0);
      for(auto &row : function_blocks[block])
        for(auto &column : row)
          for(size_t function_index(0); function_index != column.size();
              ++function_index)
            {
              auto &f(column[function_index]);
              max_delta = El::Min(max_delta, f.max_delta);
              max_degree = std::max(max_degree, f.chebyshev_coeffs.size());
              for(auto &coeff : f.chebyshev_coeffs)
                {
                  block_scale = std::max(
                    block_scale, El::Abs(coeff * weights[function_index]));
                }
            }

      const El::BigFloat block_epsilon(block_scale
                                       * El::limits::Epsilon<El::BigFloat>());

      // Preadd the coefficients
      std::vector<std::vector<Function>> summed_functions(
        function_blocks[block].size());
      for(size_t matrix_row(0); matrix_row != summed_functions.size();
          ++matrix_row)
        {
          summed_functions[matrix_row].reserve(
            function_blocks[block][matrix_row].size());
          for(size_t matrix_column(0);
              matrix_column != function_blocks[block][matrix_row].size();
              ++matrix_column)
            {
              auto &summed(summed_functions[matrix_row].emplace_back());
              // TODO: This only works if all have the same degree
              // and max_delta
              summed.max_delta = max_delta;
              summed.chebyshev_coeffs.resize(max_degree);
              auto &y_vector(
                function_blocks[block][matrix_row][matrix_column]);
              for(size_t y_index(0); y_index != y_vector.size(); ++y_index)
                {
                  for(size_t coeff_index(0);
                      coeff_index != y_vector[y_index].chebyshev_coeffs.size();
                      ++coeff_index)
                    {
                      summed.chebyshev_coeffs[coeff_index]
                        += weights[y_index]
                           * y_vector[y_index].chebyshev_coeffs[coeff_index];
                    }
                }
            }
        }

      // 1/128 should be a small enough relative error so that we are
      // in the regime of convergence.  Then the error estimates will
      // work
      Mesh mesh(
        *(points.at(block).begin()), max_delta,
        [&](const El::BigFloat &x) {
          return eval_summed(epsilon, infinity, summed_functions, x);
        },
        mesh_threshold, block_epsilon);

      new_points.at(block).clear();
      for(auto &point : get_new_points(mesh, block_epsilon))
        {
          if(points.at(block).count(point) == 0)
            {
              new_points.at(block).push_back(point);
              ++num_new_points.at(block);
            }
        }
    }
  El::mpi::AllReduce(num_new_points.data(), num_new_points.size(),
                     El::mpi::SUM, El::mpi::COMM_WORLD);

  for(size_t block(0); block != num_blocks; ++block)
    {
      new_points.at(block).resize(num_new_points.at(block));
      El::mpi::Broadcast(new_points.at(block).data(), num_new_points.at(block),
                         block % num_procs, El::mpi::COMM_WORLD);
    }
  has_new_points = (find_if(num_new_points.begin(), num_new_points.end(),
                            [](const size_t &n) { return n != 0; })
                    != num_new_points.end());
}
