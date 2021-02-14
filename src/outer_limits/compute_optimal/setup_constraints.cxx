#include "power_prefactor.hxx"
#include "poles_prefactor.hxx"
#include "../../sdp_read.hxx"

void setup_constraints(const size_t &max_index,
                       const size_t &num_blocks,
                       const El::BigFloat &infinity,
                       const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                       const std::vector<El::BigFloat> &normalization,
                       const std::vector<std::set<El::BigFloat>> &points,
                       std::vector<std::vector<El::BigFloat>> &primal_objective_c,
                       std::vector<El::Matrix<El::BigFloat>> &free_var_matrix)
{
  for(size_t block(0); block != num_blocks; ++block)
    {
      const int64_t max_degree([&]() {
        int64_t result(0);
        for(auto &row : matrices[block].polynomials)
          for(auto &column : row)
            for(auto &poly : column)
              {
                result = std::max(result, poly.degree());
              }
        return result;
      }());

      // Evalueate B and c at each point
      for(auto &x : points.at(block))
        {
          El::BigFloat prefactor([&]() {
            if(x == infinity)
              {
                return El::BigFloat(1);
              }
            return power_prefactor(matrices[block].damped_rational.base, x)
              * poles_prefactor(matrices[block].damped_rational.poles,
                                x);
          }());
          const size_t dim(matrices[block].polynomials.size());
          free_var_matrix.emplace_back(
                                       dim * (dim + 1) / 2,
                                       matrices[block].polynomials.at(0).at(0).size() - 1);
          auto &free_var(free_var_matrix.back());

          primal_objective_c.emplace_back();
          auto &primal(primal_objective_c.back());

          size_t flattened_matrix_row(0);
          for(size_t matrix_row(0); matrix_row != dim; ++matrix_row)
            for(size_t matrix_column(0); matrix_column <= matrix_row;
                ++matrix_column)
              {
                // Special case infinity by only looking at the
                // coefficients of the largest degree polynomials
                if(x == infinity)
                  {
                    if(matrices[block]
                       .polynomials.at(matrix_row)
                       .at(matrix_column)
                       .at(max_index)
                       .degree()
                       < max_degree)
                      {
                        primal.push_back(0);
                      }
                    else
                      {
                        primal.push_back(matrices[block]
                                         .polynomials.at(matrix_row)
                                         .at(matrix_column)
                                         .at(max_index)
                                         .coefficients.at(max_degree)
                                         / normalization.at(max_index));
                      }
                    auto &primal_constant(primal.back());
                    for(size_t column(0);
                        column != size_t(free_var.Width()); ++column)
                      {
                        const size_t index(column
                                           + (column < max_index ? 0 : 1));
                        if(matrices[block]
                           .polynomials.at(matrix_row)
                           .at(matrix_column)
                           .at(index)
                           .degree()
                           < max_degree)
                          {
                            free_var(flattened_matrix_row, column)
                              = primal_constant * normalization.at(index);
                          }
                        else
                          {
                            free_var(flattened_matrix_row, column)
                              = primal_constant * normalization.at(index)
                              - matrices[block]
                              .polynomials.at(matrix_row)
                              .at(matrix_column)
                              .at(index)
                              .coefficients.at(max_degree);
                          }
                      }
                  }
                else
                  {
                    primal.push_back(prefactor
                                     * matrices[block]
                                     .polynomials.at(matrix_row)
                                     .at(matrix_column)
                                     .at(max_index)(x)
                                     / normalization.at(max_index));

                    auto &primal_constant(primal.back());
                    for(size_t column(0);
                        column != size_t(free_var.Width()); ++column)
                      {
                        const size_t index(column
                                           + (column < max_index ? 0 : 1));
                        free_var(flattened_matrix_row, column)
                          = primal_constant * normalization.at(index)
                          - prefactor
                          * matrices[block]
                          .polynomials.at(matrix_row)
                          .at(matrix_column)
                          .at(index)(x);
                      }
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
