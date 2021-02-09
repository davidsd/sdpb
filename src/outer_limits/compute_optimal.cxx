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
                const SDP_Solver_Parameters &parameters_in)
{
  SDP_Solver_Parameters parameters(parameters_in);

  size_t num_weights(normalization.size());

  const size_t num_blocks(matrices.size());
  std::vector<El::BigFloat> weights(num_weights, 0);
  std::vector<std::set<El::BigFloat>> points(num_blocks);
  std::vector<std::vector<El::BigFloat>> new_points(num_blocks);

  // GMP does not have a special infinity value, so we use max double.
  const El::BigFloat zero(0), infinity(std::numeric_limits<double>::max());
  // Need to have a point at zero and infinity
  const El::BigFloat min_x(0), max_x(infinity);
  for(size_t block(0); block < num_blocks; ++block)
    {
      points.at(block).emplace(min_x);

      // points.at(block).emplace(0.1);
      // points.at(block).emplace(1);
      
      for(double x(1 / 64.0); x < 64; x *= 1.13878863476 * 1.13878863476)
        points.at(block).emplace(x);

      new_points.at(block).emplace_back(max_x);
    }

  parameters.duality_gap_threshold = 1.1;
  while(parameters.duality_gap_threshold > parameters_in.duality_gap_threshold)
    {
      if(El::mpi::Rank() == 0)
        {
          std::cout << "Threshold: " << parameters.duality_gap_threshold
                    << "\n";
        }

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

      if(El::mpi::Rank() == 0)
        {
          std::cout << "num_constraints: " << num_constraints << "\n";
        }

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

      {
        // Rescale EVERYTHING to 1
        El::BigFloat objective_scale(El::Abs(objective_const));
        for(auto &block : free_var_matrix)
          {
            for(int64_t row(0); row != block.Height(); ++row)
              for(int64_t column(0); column < block.Width(); ++column)
                {
                  block(row, column)
                    *= (objective_scale / dual_objective_b[column]);
                }
          }
        std::fill(dual_objective_b.begin(), dual_objective_b.end(),
                  El::BigFloat(1.0));
        objective_const = (objective_const > 0 ? 1 : -1);
      }

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

      if(El::mpi::Rank() == 0)
        {
          set_stream_precision(std::cout);
          std::cout << "-----" << reason << "-----\n"
                    << '\n'
                    << "primalObjective = " << solver.primal_objective << '\n'
                    << "dualObjective   = " << solver.dual_objective << '\n'
                    << "dualityGap      = " << solver.duality_gap << '\n'
                    << "primalError     = " << solver.primal_error() << '\n'
                    << "dualError       = " << solver.dual_error << '\n'
                    << '\n';
        }

      // if(reason != SDP_Solver_Terminate_Reason::PrimalDualOptimal)
      //   {
      //     std::stringstream ss;
      //     ss << "Can not find solution: " << reason;
      //     throw std::runtime_error(ss.str());
      //   }

      // y is duplicated among cores, so only need to print out copy on
      // the root node.
      // THe weight at max_index is determined by the normalization condition
      // dot(norm,weights)=1
      El::DistMatrix<El::BigFloat> y(solver.y.blocks.at(0));
      El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0), sdp.yp_to_y,
               solver.y.blocks.at(0), El::BigFloat(0.0), y);

      weights.at(max_index) = 1;
      for(size_t block_row(0); block_row != size_t(y.Height()); ++block_row)
        {
          const size_t index(block_row + (block_row < max_index ? 0 : 1));
          weights.at(index) = y.Get(block_row, 0);
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
      bool has_new_points(false);
      for(size_t block(0); block != num_blocks; ++block)
        {
          // 0.01 should be a small enough relative error so that we are
          // in the regime of convergence.  Then the error estimates will
          // work
          Mesh mesh(
            *(points.at(block).begin()), El::BigFloat(100),
            [&](const El::BigFloat &x) {
              return eval_weighted(matrices[block], x, weights);
            },
            (1.0 / 128));
          new_points.at(block) = get_new_points(mesh);
          for(auto &point : new_points.at(block))
            {
              has_new_points
                = has_new_points || (points.at(block).count(point) == 0);
            }
        }
      if(!has_new_points)
        {
          parameters.duality_gap_threshold *= (1.0 / 8);
        }
    }
  return weights;
}
