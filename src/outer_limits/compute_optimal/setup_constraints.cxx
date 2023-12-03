#include "outer_limits/Function.hxx"

#include <El.hpp>

#include <vector>

void setup_constraints(
  const size_t &max_index, const size_t &num_blocks,
  const El::BigFloat &epsilon, const El::BigFloat &infinity,
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<std::set<El::BigFloat>> &points,
  std::vector<std::vector<El::BigFloat>> &primal_objective_c,
  std::vector<El::Matrix<El::BigFloat>> &free_var_matrix)
{
  for(size_t block(0); block != num_blocks; ++block)
    {
      // Evalueate B and c at each point
      for(auto &x : points.at(block))
        {
          const size_t dim(function_blocks[block].size());
          free_var_matrix.emplace_back(
            dim * (dim + 1) / 2,
            function_blocks[block].at(0).at(0).size() - 1);
          auto &free_var(free_var_matrix.back());

          primal_objective_c.emplace_back();
          auto &primal(primal_objective_c.back());

          size_t flattened_matrix_row(0);
          for(size_t matrix_row(0); matrix_row != dim; ++matrix_row)
            for(size_t matrix_column(0); matrix_column <= matrix_row;
                ++matrix_column)
              {
                primal.push_back(function_blocks[block]
                                   .at(matrix_row)
                                   .at(matrix_column)
                                   .at(max_index)
                                   .eval(epsilon, infinity, x)
                                 / normalization.at(max_index));
                auto &primal_constant(primal.back());
                for(size_t column(0); column != size_t(free_var.Width());
                    ++column)
                  {
                    const size_t index(column + (column < max_index ? 0 : 1));
                    free_var(flattened_matrix_row, column)
                      = primal_constant * normalization.at(index)
                        - function_blocks[block]
                            .at(matrix_row)
                            .at(matrix_column)
                            .at(index)
                            .eval(epsilon, infinity, x);
                  }
                ++flattened_matrix_row;
              }
          // Rescale rows by the largest element in the row
          const El::BigFloat scaling([&]() {
            El::BigFloat max_value(0);
            size_t flattened_matrix_row(0);
            for(size_t matrix_row(0); matrix_row != dim; ++matrix_row)
              for(size_t matrix_column(0); matrix_column <= matrix_row;
                  ++matrix_column)
                {
                  max_value = std::max(
                    max_value, El::Abs(primal.at(flattened_matrix_row)));
                  for(size_t column(0); column != size_t(free_var.Width());
                      ++column)
                    {
                      max_value = std::max(
                        max_value,
                        El::Abs(free_var(flattened_matrix_row, column)));
                    }
                  ++flattened_matrix_row;
                }
            return 1 / max_value;
          }());
          {
            size_t flattened_matrix_row(0);
            for(size_t matrix_row(0); matrix_row != dim; ++matrix_row)
              for(size_t matrix_column(0); matrix_column <= matrix_row;
                  ++matrix_column)
                {
                  primal.at(flattened_matrix_row) *= scaling;
                  for(size_t column(0); column != size_t(free_var.Width());
                      ++column)
                    {
                      free_var(flattened_matrix_row, column) *= scaling;
                    }
                  ++flattened_matrix_row;
                }
          }
        }
    }
}
