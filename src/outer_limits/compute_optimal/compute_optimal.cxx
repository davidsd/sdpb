#include "../Outer_Parameters.hxx"
#include "Mesh.hxx"
#include "setup_constraints.hxx"

#include "../../ostream_set.hxx"
#include "../../ostream_map.hxx"
#include "../../ostream_vector.hxx"
#include "../../set_stream_precision.hxx"

namespace
{
  void copy_matrix(const El::Matrix<El::BigFloat> &source,
                   El::DistMatrix<El::BigFloat> &destination)
  {
    for(int64_t row(0); row < destination.LocalHeight(); ++row)
      {
        int64_t global_row(destination.GlobalRow(row));
        for(int64_t column(0); column < destination.LocalWidth(); ++column)
          {
            int64_t global_column(destination.GlobalCol(column));
            destination.SetLocal(row, column,
                                 source(global_row, global_column));
          }
      }
  }

  void
  copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
              El::Matrix<El::BigFloat> &destination)
  {
    destination.Resize(source.LocalHeight(), source.LocalWidth());
    for(int64_t row(0); row < source.LocalHeight(); ++row)
      {
        for(int64_t column(0); column < source.LocalWidth(); ++column)
          {
            destination(row, column) = source.GetLocal(row, column);
          }
      }
  }
}

void compute_y_transform(
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<std::set<El::BigFloat>> &points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const Outer_Parameters &parameters, const size_t &max_index,
  const El::Grid &global_grid,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  El::BigFloat &primal_c_scale);

El::BigFloat
eval_summed(const El::BigFloat &infinity,
            const std::vector<std::vector<Function>> &summed_functions,
            const El::BigFloat &x);

std::vector<El::BigFloat>
get_new_points(const Mesh &mesh, const El::BigFloat &block_epsilon);

void save_checkpoint(const boost::filesystem::path &checkpoint_directory,
                     const Verbosity &verbosity,
                     const boost::property_tree::ptree &parameter_properties,
                     const El::Matrix<El::BigFloat> &y,
                     const std::vector<std::set<El::BigFloat>> &points,
                     const El::BigFloat &infinity,
                     boost::optional<int64_t> &backup_generation,
                     int64_t &current_generation);

std::vector<El::BigFloat> compute_optimal(
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<std::vector<El::BigFloat>> &initial_points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const Outer_Parameters &parameters_in)
{
  if(initial_points.size() != function_blocks.size())
    {
      throw std::runtime_error(
        "Size are different: Positive_Matrix_With_Prefactor: "
        + std::to_string(function_blocks.size())
        + ", initial points: " + std::to_string(initial_points.size()));
    }
  Outer_Parameters parameters(parameters_in);

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
      points.at(block).emplace(infinity);
    }

  // TODO: This is duplicated from sdp2input/write_output/write_output.cxx
  size_t max_index(0);
  El::BigFloat max_normalization(0);
  for(size_t index(0); index != normalization.size(); ++index)
    {
      const El::BigFloat element(Abs(normalization[index]));
      if(element > max_normalization)
        {
          max_normalization = element;
          max_index = index;
        }
    }

  const El::Grid global_grid;
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> yp_to_y_star(global_grid),
    dual_objective_b_star(global_grid);
  El::BigFloat primal_c_scale;
  compute_y_transform(function_blocks, points, objectives, normalization,
                      parameters, max_index, global_grid, yp_to_y_star,
                      dual_objective_b_star, primal_c_scale);

  El::Matrix<El::BigFloat> y_saved(yp_to_y_star.Height(), 1);
  El::Zero(y_saved);

  boost::optional<int64_t> backup_generation;
  int64_t current_generation(0);
  parameters.solver.duality_gap_threshold = 1.1;
  while(parameters.solver.duality_gap_threshold
        > parameters_in.solver.duality_gap_threshold)
    {
      std::map<size_t, size_t> new_to_old;
      size_t num_constraints(0), old_index(0);
      std::vector<size_t> matrix_dimensions;
      for(size_t block(0); block != num_blocks; ++block)
        {
          for(size_t offset(0); offset != points.at(block).size(); ++offset)
            {
              new_to_old.emplace(num_constraints + offset, old_index + offset);
            }
          old_index += points.at(block).size();
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
      if(rank == 0 && parameters.verbosity >= Verbosity::regular)
        {
          std::cout << "num_constraints: " << num_constraints << "\n";
        }

      std::vector<std::vector<El::BigFloat>> primal_objective_c;
      primal_objective_c.reserve(num_constraints);
      std::vector<El::Matrix<El::BigFloat>> free_var_matrix;
      free_var_matrix.reserve(num_constraints);

      setup_constraints(max_index, num_blocks, infinity, function_blocks,
                        normalization, points, primal_objective_c,
                        free_var_matrix);

      const El::BigFloat objective_const(objectives.at(max_index)
                                         / normalization.at(max_index));

      Block_Info block_info(matrix_dimensions, parameters.verbosity);

      El::Grid grid(block_info.mpi_comm.value);

      SDP sdp(objective_const, primal_objective_c, free_var_matrix,
              yp_to_y_star, dual_objective_b_star, primal_c_scale, block_info,
              grid);

      SDP_Solver solver(parameters.solver, parameters.verbosity,
                        parameters.require_initial_checkpoint, block_info,
                        grid, sdp.dual_objective_b.Height());

      for(auto &y_block : solver.y.blocks)
        {
          copy_matrix(y_saved, y_block);
        }

      boost::property_tree::ptree parameter_properties(
        to_property_tree(parameters));
      bool has_new_points(false);
      while(!has_new_points
            && parameters.solver.duality_gap_threshold
                 > parameters_in.solver.duality_gap_threshold)
        {
          if(rank == 0 && parameters.verbosity >= Verbosity::regular)
            {
              std::cout << "Threshold: "
                        << parameters.solver.duality_gap_threshold << "\n";
            }

          Timers timers(parameters.verbosity >= Verbosity::debug);
          SDP_Solver_Terminate_Reason reason
            = solver.run(parameters.solver, parameters.verbosity,
                         parameter_properties, block_info, sdp, grid, timers);

          if(rank == 0 && parameters.verbosity >= Verbosity::regular)
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

          if(reason == SDP_Solver_Terminate_Reason::MaxComplementarityExceeded
             || reason == SDP_Solver_Terminate_Reason::MaxIterationsExceeded
             || reason == SDP_Solver_Terminate_Reason::MaxRuntimeExceeded)
            {
              std::stringstream ss;
              ss << "Can not find solution: " << reason;
              throw std::runtime_error(ss.str());
            }

          // y is duplicated among cores, so only need to print out copy on
          // the root node.
          // THe weight at max_index is determined by the normalization
          // condition dot(norm,weights)=1
          El::DistMatrix<El::BigFloat> yp(dual_objective_b_star.Height(), 1,
                                          yp_to_y_star.Grid());
          El::Zero(yp);
          El::DistMatrix<El::BigFloat> y(yp);
          El::DistMatrix<El::BigFloat, El::STAR, El::STAR> yp_star(
            solver.y.blocks.at(0));
          for(int64_t row(0); row != yp.LocalHeight(); ++row)
            {
              int64_t global_row(yp.GlobalRow(row));
              for(int64_t column(0); column != yp.LocalWidth(); ++column)
                {
                  int64_t global_column(yp.GlobalCol(column));
                  yp.SetLocal(row, column,
                              yp_star.GetLocal(global_row, global_column));
                }
            }
          El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0), yp_to_y_star,
                   yp, El::BigFloat(0.0), y);
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
          if(rank == 0 && parameters.verbosity >= Verbosity::regular)
            {
              set_stream_precision(std::cout);
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
              // TODO: These can both be precomputed
              El::BigFloat max_delta(infinity), block_scale(0);
              size_t max_degree(0);
              for(auto &row : function_blocks[block])
                for(auto &column : row)
                  for(size_t function_index(0);
                      function_index != column.size(); ++function_index)
                    {
                      auto &f(column[function_index]);
                      max_delta = El::Min(max_delta, f.max_delta);
                      max_degree
                        = std::max(max_degree, f.chebyshev_coeffs.size());
                      for(auto &coeff : f.chebyshev_coeffs)
                        {
                          block_scale = std::max(
                            block_scale,
                            El::Abs(coeff * weights[function_index]));
                        }
                    }

              const El::BigFloat block_epsilon(
                block_scale * El::limits::Epsilon<El::BigFloat>());

              // Preadd the coefficients
              std::vector<std::vector<Function>> summed_functions(
                function_blocks[block].size());
              El::BigFloat zero(0);
              for(size_t matrix_row(0); matrix_row != summed_functions.size();
                  ++matrix_row)
                {
                  summed_functions[matrix_row].reserve(
                    function_blocks[block][matrix_row].size());
                  for(size_t matrix_column(0);
                      matrix_column
                      != function_blocks[block][matrix_row].size();
                      ++matrix_column)
                    {
                      auto &summed(
                        summed_functions[matrix_row].emplace_back());
                      // TODO: This only works if all have the same degree
                      // and max_delta
                      summed.max_delta = max_delta;
                      summed.infinity_value = zero;
                      summed.chebyshev_coeffs.resize(max_degree);
                      auto &y_vector(
                        function_blocks[block][matrix_row][matrix_column]);
                      for(size_t y_index(0); y_index != y_vector.size();
                          ++y_index)
                        {
                          for(size_t coeff_index(0);
                              coeff_index
                              != y_vector[y_index].chebyshev_coeffs.size();
                              ++coeff_index)
                            {
                              summed.chebyshev_coeffs[coeff_index]
                                += weights[y_index]
                                   * y_vector[y_index]
                                       .chebyshev_coeffs[coeff_index];
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
                  return eval_summed(infinity, summed_functions, x);
                },
                (1.0 / 128), block_epsilon);

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
              parameters.solver.duality_gap_threshold *= (1.0 / 8);
            }
        }
      El::DistMatrix<El::BigFloat, El::STAR, El::STAR> y_star(
        solver.y.blocks.front());
      copy_matrix(y_star, y_saved);
      save_checkpoint(parameters.solver.checkpoint_out, parameters.verbosity,
                      parameter_properties, y_saved, points, infinity,
                      backup_generation, current_generation);
    }
  return weights;
}
