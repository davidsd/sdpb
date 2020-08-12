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

  std::cout << "blocks: " << num_blocks << "\n";
  const El::BigFloat min_x(0), max_x(10);
  for(size_t block(0); block < num_blocks; ++block)
    {
      points.at(block).emplace(min_x);
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
          for(auto &point: points.at(block))
            {
              matrix_dimensions.push_back(matrices[block].polynomials.size());
            }
          std::cout << "points: " << block << " " << points.at(block) << "\n";
        }


      std::cout << "num_constraints: " << num_constraints << "\n";
      std::cout << "matrix_dimensions: " << matrix_dimensions << "\n";

      Block_Info block_info(matrix_dimensions, parameters.procs_per_node,
                            parameters.proc_granularity, parameters.verbosity);
      El::Grid grid(block_info.mpi_comm.value);
      std::vector<El::BigFloat> prefactors;
      prefactors.reserve(num_constraints);
      std::vector<std::vector<El::BigFloat>> primal_objective_c;
      primal_objective_c.reserve(num_constraints);
      std::vector<El::Matrix<El::BigFloat>> free_var_matrix;
      free_var_matrix.reserve(num_constraints);
      for(size_t block(0); block != num_blocks; ++block)
        {
          for(auto &x : points.at(block))
            {
              prefactors.push_back(
                power_prefactor(matrices[block].damped_rational.base, x)
                * poles_prefactor(matrices[block].damped_rational.poles, x));
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
                    primal.push_back(prefactor
                                     * matrices[block]
                                         .polynomials.at(matrix_row)
                                         .at(matrix_column)
                                         .at(0)(x));

                    for(int64_t column(0); column != free_var.Width();
                        ++column)
                      {
                        free_var(flattened_matrix_row, column)
                          = -prefactor
                            * matrices[block]
                                .polynomials.at(matrix_row)
                                .at(matrix_column)
                                .at(column + 1)(x);
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

      // y is duplicated among cores, so only need to print out copy on
      // the root node.
      // THe first weight is implicitly 1.
      for(int64_t block_row(0); block_row != solver.y.blocks.at(0).Height();
          ++block_row)
        {
          weights.at(block_row + 1) = solver.y.blocks.at(0).Get(block_row, 0);
        }

      std::cout << "weight: " << weights << "\n";
      for(size_t block(0); block != num_blocks; ++block)
        {
          // 0.01 should be a small enough relative error so that we are
          // in the regime of convergence.  Then the error estimates will
          // work
          Mesh mesh(*(points.at(block).begin()), *(points.at(block).rbegin()),
                    [&](const El::BigFloat &x) {
                      return eval_weighted(matrices[block], x, weights);
                    },
                    0.01);
          new_points.at(block) = get_new_points(mesh);
          // std::cout << "new: " << block << " " << new_points.at(block) <<
          // "\n";
          has_new_points = has_new_points || !new_points.at(block).empty();
        }
    }
  // std::cout.precision(precision / 3.3);
  std::cout << "weights: " << weights << "\n";
  return weights;
}
