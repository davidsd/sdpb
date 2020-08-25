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
  std::vector<El::BigFloat> weights(num_weights, 1);
  std::vector<std::set<El::BigFloat>> points(num_blocks);
  std::vector<std::vector<El::BigFloat>> new_points(num_blocks);

  // Need to have a point at zero and infinity
  // GMP does not have a special infinity value, so we use max double.
  const El::BigFloat infinity(std::numeric_limits<double>::max());
  const El::BigFloat min_x(0), max_x(infinity);
  for(size_t block(0); block < num_blocks; ++block)
    {
      points.at(block).emplace(min_x);
      points.at(block).emplace(1);
      // for(double x(min_x); x < max_x; x *= 4)
      //   points.at(block).emplace(x);

      new_points.at(block).emplace_back(max_x);
    }

  bool has_new_points(true);

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

      if(El::mpi::Rank() == 0)
        {
          std::cout << "num_constraints: " << num_constraints << "\n";
        }
      // std::cout << "matrix_dimensions: " << matrix_dimensions << "\n";

      Block_Info block_info(matrix_dimensions, parameters.procs_per_node,
                            parameters.proc_granularity, parameters.verbosity);
      El::Grid grid(block_info.mpi_comm.value);
      std::vector<El::BigFloat> prefactors;
      prefactors.reserve(num_constraints);
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
      int64_t max_index(
        std::distance(normalization.begin(), max_normalization));

      for(size_t block(0); block != num_blocks; ++block)
        {
          for(auto &x : points.at(block))
            {
              if(x == infinity)
                {
                  prefactors.push_back(1);
                }
              else
                {
                  prefactors.push_back(
                    power_prefactor(matrices[block].damped_rational.base, x)
                    * poles_prefactor(matrices[block].damped_rational.poles,
                                      x));
                }
              auto &prefactor(prefactors.back());
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
                        int64_t max_degree(0);
                        for(auto &poly : matrices[block]
                                           .polynomials.at(matrix_row)
                                           .at(matrix_column))
                          max_degree = std::max(max_degree, poly.degree());
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
                        for(int64_t column(0); column != free_var.Width();
                            ++column)
                          {
                            const int64_t index(
                              column + (column < max_index ? 0 : 1));
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
                        for(int64_t column(0); column != free_var.Width();
                            ++column)
                          {
                            const int64_t index(
                              column + (column < max_index ? 0 : 1));
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
            }
        }

      SDP sdp(objectives, normalization, prefactors, primal_objective_c,
              free_var_matrix, block_info, grid);

      SDP_Solver solver(parameters, block_info, grid,
                        sdp.dual_objective_b.Height());

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
      for(int64_t block_row(0); block_row != solver.y.blocks.at(0).Height();
          ++block_row)
        {
          const int64_t index(block_row + (block_row < max_index ? 0 : 1));
          weights.at(index) = solver.y.blocks.at(0).Get(block_row, 0);
          weights.at(max_index) -= weights.at(index) * normalization.at(index);
        }
      weights.at(max_index) /= normalization.at(max_index);

      if(El::mpi::Rank() == 0)
        {
          std::cout.precision(10);
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
