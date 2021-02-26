#include "Mesh.hxx"

#include "../Function.hxx"
#include "../../sdp_solve.hxx"
#include "../../ostream_set.hxx"
#include "../../ostream_vector.hxx"
#include "../../set_stream_precision.hxx"

void setup_constraints_functions(
  const size_t &max_index, const size_t &num_blocks,
  const El::BigFloat &infinity,
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<std::set<El::BigFloat>> &points,
  std::vector<std::vector<El::BigFloat>> &primal_objective_c,
  std::vector<El::Matrix<El::BigFloat>> &free_var_matrix);

El::BigFloat eval_weighted_functions(
  const El::BigFloat &infinity,
  const std::vector<std::vector<std::vector<Function>>> &function_blocks,
  const El::BigFloat &x, const std::vector<El::BigFloat> &weights);

std::vector<El::BigFloat> get_new_points(const Mesh &mesh);

std::vector<El::BigFloat> compute_optimal_functions(
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<std::vector<El::BigFloat>> &initial_points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const SDP_Solver_Parameters &parameters_in)
{
  if(initial_points.size() != function_blocks.size())
    {
      throw std::runtime_error(
        "Size are different: Positive_Matrix_With_Prefactor: "
        + std::to_string(function_blocks.size())
        + ", initial points: " + std::to_string(initial_points.size()));
    }
  SDP_Solver_Parameters parameters(parameters_in);

  const size_t rank(El::mpi::Rank()), num_procs(El::mpi::Size()),
    num_weights(normalization.size());

  const size_t num_blocks(initial_points.size());
  std::vector<El::BigFloat> weights(num_weights, 0);
  std::vector<std::set<El::BigFloat>> points(num_blocks);
  std::vector<std::vector<El::BigFloat>> new_points(num_blocks);

  // GMP does not have a special infinity value, so we use max double.
  const El::BigFloat infinity(std::numeric_limits<double>::max());
  // Use the input points and add inifinty
  for(size_t block(0); block < num_blocks; ++block)
    {
      for(auto &point : initial_points.at(block))
        {
          points.at(block).emplace(point);
        }
      new_points.at(block).emplace_back(infinity);
    }

  parameters.duality_gap_threshold = 1.1;
  while(parameters.duality_gap_threshold > parameters_in.duality_gap_threshold)
    {
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
                                   function_blocks[block].size());
          if(rank == 0 && parameters.verbosity >= Verbosity::debug)
            {
              std::cout << "points: " << block << " " << points.at(block)
                        << "\n";
            }
        }

      if(rank == 0)
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

      setup_constraints_functions(max_index, num_blocks, infinity,
                                  function_blocks, normalization, points,
                                  primal_objective_c, free_var_matrix);

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

      bool has_new_points(false);
      while(!has_new_points
            && parameters.duality_gap_threshold
                 > parameters_in.duality_gap_threshold)
        {
          if(rank == 0)
            {
              std::cout << "Threshold: " << parameters.duality_gap_threshold
                        << "\n";
            }

          Timers timers(parameters.verbosity >= Verbosity::debug);
          SDP_Solver_Terminate_Reason reason
            = solver.run(parameters, block_info, sdp, grid, timers);

          if(rank == 0)
            {
              set_stream_precision(std::cout);
              std::cout << "-----" << reason << "-----\n"
                        << '\n'
                        << "primalObjective = " << solver.primal_objective
                        << '\n'
                        << "dualObjective   = " << solver.dual_objective
                        << '\n'
                        << "dualityGap      = " << solver.duality_gap << '\n'
                        << "primalError     = " << solver.primal_error()
                        << '\n'
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
          // THe weight at max_index is determined by the normalization
          // condition dot(norm,weights)=1
          El::DistMatrix<El::BigFloat> y(sdp.yp_to_y.Grid());
          y.Resize(dual_objective_b.size(), 1);
          El::Zero(y);
          El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0), sdp.yp_to_y,
                   solver.y.blocks.at(0), El::BigFloat(0.0), y);
          El::DistMatrix<El::BigFloat, El::STAR, El::STAR> y_star(y);

          weights.at(max_index) = 1;
          for(size_t block_row(0); block_row != size_t(y_star.LocalHeight());
              ++block_row)
            {
              const size_t index(block_row + (block_row < max_index ? 0 : 1));
              weights.at(index) = y_star.GetLocalCRef(block_row, 0);
              weights.at(max_index)
                -= weights.at(index) * normalization.at(index);
            }
          weights.at(max_index) /= normalization.at(max_index);
          if(rank == 0)
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
          std::vector<size_t> num_new_points(num_blocks, 0);
          for(size_t block(rank); block < num_blocks; block += num_procs)
            {
              El::BigFloat max_delta(infinity);
              for(auto &row: function_blocks[block])
                for(auto &column: row)
                  for(auto &f: column)
                    {
                      max_delta=El::Min(max_delta,f.max_delta);
                    }
                
              // 1/128 should be a small enough relative error so that we are
              // in the regime of convergence.  Then the error estimates will
              // work

              Mesh mesh(
                *(points.at(block).begin()), max_delta,
                [&](const El::BigFloat &x) {
                  return eval_weighted_functions(
                    infinity, function_blocks[block], x, weights);
                },
                (1.0 / 128));
              std::vector<El::BigFloat> candidates(get_new_points(mesh));
              new_points.at(block).clear();
              for(auto &point : candidates)
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
              El::mpi::Broadcast(new_points.at(block).data(),
                                 num_new_points.at(block), block % num_procs,
                                 El::mpi::COMM_WORLD);
            }
          has_new_points
            = (find_if(num_new_points.begin(), num_new_points.end(),
                       [](const size_t &n) { return n != 0; })
               != num_new_points.end());
          if(!has_new_points)
            {
              parameters.duality_gap_threshold *= (1.0 / 8);
            }
        }
    }
  return weights;
}
