#include "get_new_points.hxx"

#include "eval_weighted.hxx"
#include "poles_prefactor.hxx"
#include "power_prefactor.hxx"

#include "../sdp_solve.hxx"
#include "../ostream_set.hxx"

std::vector<El::BigFloat>
compute_optimal(const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                const std::vector<El::BigFloat> &objectives,
                const std::vector<El::BigFloat> &normalization,
                const SDP_Solver_Parameters &parameters)
{
  size_t num_weights(normalization.size());

  const size_t num_blocks(matrices.size());
  std::vector<El::BigFloat> weights(num_weights, 0);
  std::vector<std::set<El::BigFloat>> points(num_blocks);
  std::vector<std::vector<El::BigFloat>> new_points(num_blocks);

  // Need to have a point at zero and infinity
  // GMP does not have a special infinity value, so we use max double.
  const El::BigFloat infinity(std::numeric_limits<double>::max());
  const El::BigFloat min_x(0), max_x(infinity);
  for(size_t block(0); block < num_blocks; ++block)
    {
      points.at(block).emplace(min_x);

      // points.at(block).emplace(0.1);
      // points.at(block).emplace(1);

      // const int64_t num_points(64);
      // const double dx(1.0/num_points);
      // for(double x(dx); x < 1; x +=dx)
      //   points.at(block).emplace(-log(1-x));

      // const int64_t num_points(64);
      // const double dx(64.0/num_points);
      // for(double x(dx); x <= 64; x +=dx)
      //   points.at(block).emplace(x);

      for(double x(1 / 64.0); x < 64; x *= 1.13878863476 * 1.13878863476)
        points.at(block).emplace(x);

      new_points.at(block).emplace_back(max_x);
    }

  bool has_new_points(true);

  bool is_first_iteration(true);
  while(has_new_points)
    {
      has_new_points = false;
      size_t num_constraints(0);
      std::vector<size_t> matrix_dimensions;
      for(size_t block(0); block != num_blocks; ++block)
        {
          for(auto &point : new_points.at(block))
            {
              points.at(block).emplace(point);
            }
          num_constraints += points.at(block).size();
          matrix_dimensions.insert(matrix_dimensions.end(),
                                   points.at(block).size(),
                                   matrices[block].polynomials.size());
          if(El::mpi::Rank() == 0)
            {
              std::cout << "points: " << block << " " << points.at(block)
                        << "\n";
            }
        }

      // if(El::mpi::Rank() == 0)
      //   {
      //     std::cout << "num_constraints: " << num_constraints << "\n";
      //   }

      std::vector<std::vector<El::BigFloat>> primal_objective_c;
      primal_objective_c.reserve(num_constraints);
      std::vector<El::Matrix<El::BigFloat>> free_var_matrix;
      free_var_matrix.reserve(num_constraints);

      // TODO: This is duplicated from sdp2input/write_output/write_output.cxx
      auto max_normalization(normalization.begin());
      for(auto n(normalization.begin()); n != normalization.end(); ++n)
        {
          if(Abs(*n) > Abs(*max_normalization))
            {
              max_normalization = n;
            }
        }
      const size_t max_index(
        std::distance(normalization.begin(), max_normalization));

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
              // Rescale rows
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

      El::BigFloat objective_const(objectives.at(max_index)
                                   / normalization.at(max_index));
      std::vector<El::BigFloat> dual_objective_b;
      dual_objective_b.reserve(normalization.size() - 1);
      for(size_t index = 0; index < normalization.size(); ++index)
        {
          if(index != max_index)
            {
              dual_objective_b.push_back(objectives.at(index)
                                         - normalization.at(index)
                                             * objective_const);
            }
        }

      // Rescale columns
      std::vector<El::BigFloat> rescaling(num_weights - 1, 0);
      for(auto &matrix : free_var_matrix)
        {
          for(size_t column(0); column < size_t(matrix.Width()); ++column)
            {
              for(size_t row(0); row < size_t(matrix.Height()); ++row)
                {
                  rescaling.at(column) = std::max(
                    rescaling.at(column), El::Abs(matrix(row, column)));
                }
            }
        }

      for(size_t index(0); index != rescaling.size(); ++index)
        {
          rescaling[index] = 1 / rescaling[index];
          dual_objective_b[index] *= rescaling[index];
        }

      for(auto &matrix : free_var_matrix)
        {
          for(size_t row(0); row != size_t(matrix.Height()); ++row)
            {
              for(size_t column(0); column != size_t(matrix.Width()); ++column)
                {
                  matrix(row, column) *= rescaling[column];
                }
            }
        }

      if(is_first_iteration)
        {
          std::vector<size_t> elements_to_keep(1, 0);
          const El::BigFloat tolerance(0.1);
          for(size_t matrix_index(0);
              matrix_index + 1 < free_var_matrix.size();)
            {
              auto &matrix(free_var_matrix[matrix_index]);
              auto &c(primal_objective_c[matrix_index]);

              size_t offset_index(matrix_index + 1);
              for(; offset_index != free_var_matrix.size(); ++offset_index)
                {
                  if(free_var_matrix[offset_index].Height() != matrix.Height())
                    {
                      elements_to_keep.push_back(offset_index);
                      break;
                    }
                  bool should_erase(true);
                  for(size_t row(0);
                      should_erase && row != size_t(matrix.Height()); ++row)
                    {
                      should_erase
                        = should_erase
                          && (El::Abs(c[row]
                                      - primal_objective_c[offset_index][row])
                              < tolerance);
                      for(size_t column(0);
                          should_erase && column != size_t(matrix.Width());
                          ++column)
                        {
                          should_erase
                            = should_erase
                              && (El::Abs(matrix(row, column)
                                          - free_var_matrix[offset_index](
                                              row, column))
                                  < tolerance);
                        }
                    }
                  if(!should_erase)
                    {
                      elements_to_keep.push_back(offset_index);
                      break;
                    }
                }
              matrix_index = offset_index;
            }

          std::vector<std::vector<El::BigFloat>> temp_c;
          temp_c.reserve(elements_to_keep.size());
          std::vector<El::Matrix<El::BigFloat>> temp_matrix;
          temp_matrix.reserve(elements_to_keep.size());
          std::vector<size_t> temp_dimensions;
          temp_dimensions.reserve(elements_to_keep.size());

          for(auto element : elements_to_keep)
            {
              temp_c.push_back(primal_objective_c[element]);
              temp_matrix.push_back(free_var_matrix[element]);
              temp_dimensions.push_back(matrix_dimensions[element]);
            }

          std::vector<std::set<El::BigFloat>> temp_points(num_blocks);
          size_t index(0), temp_num_constraints(0);
          for(size_t block(0); block != num_blocks; ++block)
            {
              for(auto &point : points.at(block))
                {
                  if(std::find(elements_to_keep.begin(),
                               elements_to_keep.end(), index)
                     != elements_to_keep.end())
                    temp_points.at(block).emplace(point);
                  ++index;
                }
              temp_num_constraints += temp_points.at(block).size();
              if(El::mpi::Rank() == 0)
                {
                  std::cout << "points: " << block << " "
                            << temp_points.at(block) << "\n";
                }
            }

          std::swap(temp_c, primal_objective_c);
          std::swap(temp_matrix, free_var_matrix);
          std::swap(temp_dimensions, matrix_dimensions);
          std::swap(temp_points, points);

          is_first_iteration = false;
          // std::cout << "dimensions: " << temp_dimensions.size() << " "
          //           << matrix_dimensions.size() << " "
          //           << temp_num_constraints << "\n";
          // std::cout << matrix_dimensions << "\n";
          // std::cout << points << "\n";
        }
      std::cout << "num_constraints: " << free_var_matrix.size() << "\n";

      Block_Info block_info(matrix_dimensions, parameters.procs_per_node,
                            parameters.proc_granularity, parameters.verbosity);
      El::Grid grid(block_info.mpi_comm.value);

      SDP sdp(objective_const, dual_objective_b, primal_objective_c,
              free_var_matrix, block_info, grid);

      SDP_Solver solver(parameters, block_info, grid,
                        sdp.dual_objective_b.Height());

      for(auto &block : solver.y.blocks)
        {
          if(block.GlobalCol(0) == 0)
            {
              for(size_t row(0); row != size_t(block.LocalHeight()); ++row)
                {
                  size_t global_row(block.GlobalRow(row));
                  const size_t index(global_row
                                     + (global_row < max_index ? 0 : 1));
                  block.SetLocal(row, 0, weights.at(index));
                }
            }
        }

      Timers timers(parameters.verbosity >= Verbosity::debug);
      SDP_Solver_Terminate_Reason reason
        = solver.run(parameters, block_info, sdp, grid, timers);

      if(reason != SDP_Solver_Terminate_Reason::PrimalDualOptimal)
        {
          std::stringstream ss;
          ss << "Can not find solution: " << reason;
          throw std::runtime_error(ss.str());
        }
      // y is duplicated among cores, so only need to print out copy on
      // the root node.
      // THe weight at max_index is determined by the normalization condition
      // dot(norm,weights)=1
      weights.at(max_index) = 1;
      for(size_t block_row(0);
          block_row != size_t(solver.y.blocks.at(0).Height()); ++block_row)
        {
          const size_t index(block_row + (block_row < max_index ? 0 : 1));
          weights.at(index)
            = solver.y.blocks.at(0).Get(block_row, 0) * rescaling[block_row];
          weights.at(max_index) -= weights.at(index) * normalization.at(index);
        }
      weights.at(max_index) /= normalization.at(max_index);

      if(El::mpi::Rank() == 0)
        {
          std::cout.precision(20);
          std::cout << "weight: " << weights << "\n";

          El::BigFloat optimal(0);
          for(size_t index(0); index < objectives.size(); ++index)
            {
              optimal += objectives[index] * weights[index];
            }
          std::cout << "optimal: " << optimal << "\n";
        }
      for(size_t block(0); block != num_blocks; ++block)
        {
          // 0.01 should be a small enough relative error so that we are
          // in the regime of convergence.  Then the error estimates will
          // work
          // Mesh mesh(*(points.at(block).begin()), *(points.at(block).rbegin()),
          Mesh mesh(*(points.at(block).begin()), El::BigFloat(100),
                    [&](const El::BigFloat &x) {
                      return eval_weighted(matrices[block], x, weights);
                    },
                    0.01);
          new_points.at(block) = get_new_points(mesh);
          for(auto &point : new_points.at(block))
            {
              has_new_points
                = has_new_points || (points.at(block).count(point) == 0);
            }
        }
    }
  // if(El::mpi::Rank() == 0)
  //   {
  //     std::cout << "weights: " << weights << "\n";
  //   }
  return weights;
}
