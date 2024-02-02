#include "setup_constraints.hxx"
#include "setup_constraints.hxx"
#include "outer_limits/Outer_Parameters.hxx"
#include "pmp/max_normalization_index.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/copy_matrix.hxx"
#include "sdpb_util/fill_weights.hxx"
#include "sdpb_util/ostream/ostream_map.hxx"
#include "sdpb_util/ostream/ostream_set.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

namespace fs = std::filesystem;

void compute_y_transform(
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<std::set<El::BigFloat>> &points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization, const Environment &env,
  const Outer_Parameters &parameters, const size_t &max_index,
  const El::Grid &global_grid,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  El::BigFloat &primal_c_scale);

boost::optional<int64_t> load_checkpoint(
  const fs::path &checkpoint_directory, const Verbosity &verbosity,
  boost::optional<int64_t> &backup_generation, int64_t &current_generation,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y_star,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  El::Matrix<El::BigFloat> &y, std::vector<std::set<El::BigFloat>> &points,
  El::BigFloat &threshold, El::BigFloat &primal_c_scale);

void find_new_points(
  const size_t &num_blocks, const size_t &rank, const size_t &num_procs,
  const El::BigFloat &mesh_threshold, const El::BigFloat &epsilon,
  const El::BigFloat &infinity,
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<El::BigFloat> &weights,
  const std::vector<std::set<El::BigFloat>> &points,
  std::vector<std::vector<El::BigFloat>> &new_points, bool &has_new_points);

void save_checkpoint(
  const fs::path &checkpoint_directory, const Verbosity &verbosity,
  const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y_star,
  const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  const El::Matrix<El::BigFloat> &y,
  const std::vector<std::set<El::BigFloat>> &points,
  const El::BigFloat &infinity, const El::BigFloat &threshold,
  const El::BigFloat &primal_c_scale,
  boost::optional<int64_t> &backup_generation, int64_t &current_generation);

std::vector<El::BigFloat> compute_optimal(
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<std::vector<El::BigFloat>> &initial_points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization, const Environment &env,
  const Outer_Parameters &parameters_in,
  const std::chrono::time_point<std::chrono::high_resolution_clock> &start_time)
{
  ASSERT(initial_points.size() == function_blocks.size(),
         "Size are different: Polynomial_Vector_Matrix: "
           + std::to_string(function_blocks.size())
           + ", initial points: " + std::to_string(initial_points.size()));
  Outer_Parameters parameters(parameters_in);
  const size_t rank(El::mpi::Rank()), num_procs(El::mpi::Size()),
    num_weights(normalization.size());

  const size_t num_blocks(initial_points.size());
  std::vector<El::BigFloat> weights(num_weights, 0);
  std::vector<std::set<El::BigFloat>> points(num_blocks);
  std::vector<std::vector<El::BigFloat>> new_points(num_blocks);

  // GMP does not have a special infinity value, so we use max double.
  const El::BigFloat infinity(std::numeric_limits<double>::max()),
    epsilon(El::limits::Epsilon<El::BigFloat>());
  // Use the input points and add inifinty
  for(size_t block(0); block < num_blocks; ++block)
    {
      points.at(block).emplace(epsilon);
      for(auto &point : initial_points.at(block))
        {
          points.at(block).emplace(point);
        }
      points.at(block).emplace(infinity);
    }

  const size_t max_index(max_normalization_index(normalization));

  const El::Grid global_grid;
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> yp_to_y_star(
    normalization.size() - 1, normalization.size() - 1, global_grid),
    dual_objective_b_star(global_grid);
  El::BigFloat primal_c_scale;

  parameters.solver.duality_gap_threshold = 1.1;
  El::Matrix<El::BigFloat> yp_saved(yp_to_y_star.Height(), 1);
  El::Zero(yp_saved);

  int64_t current_generation(0);
  boost::optional<int64_t> backup_generation;
  load_checkpoint(parameters.solver.checkpoint_in, parameters_in.verbosity,
                  backup_generation, current_generation, yp_to_y_star,
                  dual_objective_b_star, yp_saved, points,
                  parameters.solver.duality_gap_threshold, primal_c_scale);
  if(backup_generation)
    {
      El::Matrix<El::BigFloat> y(yp_saved.Height(), 1);

      El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
               yp_to_y_star.LockedMatrix(), yp_saved, El::BigFloat(0.0), y);
      fill_weights(y, max_index, normalization, weights);
    }
  else
    {
      compute_y_transform(function_blocks, points, objectives, normalization,
                          env, parameters, max_index, global_grid,
                          yp_to_y_star, dual_objective_b_star, primal_c_scale);
    }

  while(parameters.solver.duality_gap_threshold
        >= parameters_in.solver.duality_gap_threshold)
    {
      std::map<size_t, size_t> new_to_old;
      size_t num_constraints(0), old_index(0);
      std::vector<size_t> matrix_dimensions;
      for(size_t block(0); block != num_blocks; ++block)
        {
          std::stringstream ss;
          set_stream_precision(ss);
          ss << "new_points: "
             << " " << block;

          for(size_t offset(0); offset != points.at(block).size(); ++offset)
            {
              new_to_old.emplace(num_constraints + offset, old_index + offset);
            }
          old_index += points.at(block).size();
          for(auto &point : new_points.at(block))
            {
              points.at(block).emplace(point);
              ss << " " << point;
            }
          ss << "\n";
          if(El::mpi::Rank() == 0 && !new_points.at(block).empty()
             && parameters_in.verbosity >= Verbosity::regular)
            {
              std::cout << ss.str() << std::flush;
            }
          num_constraints += points.at(block).size();
          matrix_dimensions.insert(matrix_dimensions.end(),
                                   points.at(block).size(),
                                   function_blocks[block].size());
          if(rank == 0 && parameters_in.verbosity >= Verbosity::debug)
            {
              std::cout << "points: " << block << " " << points.at(block)
                        << "\n";
            }
        }
      if(rank == 0 && parameters_in.verbosity >= Verbosity::regular)
        {
          std::cout << "num_constraints: " << num_constraints << "\n";
        }

      std::vector<std::vector<El::BigFloat>> primal_objective_c;
      primal_objective_c.reserve(num_constraints);
      std::vector<El::Matrix<El::BigFloat>> free_var_matrix;
      free_var_matrix.reserve(num_constraints);

      setup_constraints(max_index, num_blocks, epsilon, infinity,
                        function_blocks, normalization, points,
                        primal_objective_c, free_var_matrix);

      const El::BigFloat objective_const(objectives.at(max_index)
                                         / normalization.at(max_index));

      Block_Info block_info(env, matrix_dimensions, parameters.verbosity);

      El::Grid grid(block_info.mpi_comm.value);

      SDP sdp(objective_const, primal_objective_c, free_var_matrix,
              yp_to_y_star, dual_objective_b_star, normalization, primal_c_scale, block_info,
              grid);

      SDP_Solver solver(parameters.solver, parameters.verbosity,
                        parameters.require_initial_checkpoint, block_info,
                        grid, sdp.dual_objective_b.Height());

      for(auto &y_block : solver.y.blocks)
        {
          copy_matrix(yp_saved, y_block);
        }

      boost::property_tree::ptree parameter_properties(
        to_property_tree(parameters));
      bool has_new_points(false);
      while(!has_new_points
            && parameters.solver.duality_gap_threshold
                 >= parameters_in.solver.duality_gap_threshold)
        {
          if(rank == 0 && parameters_in.verbosity >= Verbosity::regular)
            {
              std::cout << "Threshold: "
                        << parameters.solver.duality_gap_threshold << "\n";
            }

          Timers timers(env, parameters.verbosity >= Verbosity::debug);
          SDP_Solver_Terminate_Reason reason = solver.run(
            parameters.solver, parameters.verbosity, parameter_properties,
            block_info, sdp, grid, start_time, timers);

          for(size_t index(0); index < block_info.block_indices.size();
              ++index)
            {
              size_t block_number(0),
                point_index(block_info.block_indices[index]);
              while(point_index > points.at(block_number).size())
                {
                  point_index -= points.at(block_number).size();
                  ++block_number;
                }
              if(block_number == 2)
                {
                  auto &block(solver.Y.blocks[index * 2]);
                  if(block.Height() > 0)
                    {
                      El::BigFloat determinant(-1);
                      switch(block.Height())
                        {
                        case 1: determinant = block.Get(0, 0); break;
                        case 2:
                          determinant = block.Get(0, 0) * block.Get(1, 1)
                                        - block.Get(0, 1) * block.Get(0, 1);
                          break;
                        default:
                          RUNTIME_ERROR("too big: ", block.Height());
                          break;
                        }
                      if(determinant < 1e-16)
                        {
                          std::cout.precision(20);
                          std::cout
                            << "Y: " << El::mpi::Rank() << " " << index << " "
                            << block.Height() << " " << block_number << " "
                            << *(std::next(points.at(block_number).begin(),
                                           point_index))
                            << " " << determinant << "\n"
                            << std::flush;
                        }
                    }
                }
            }

          if(rank == 0 && parameters_in.verbosity >= Verbosity::debug)
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
             || reason == SDP_Solver_Terminate_Reason::MaxRuntimeExceeded
             || reason == SDP_Solver_Terminate_Reason::PrimalStepTooSmall
             || reason == SDP_Solver_Terminate_Reason::DualStepTooSmall)
            {
              RUNTIME_ERROR("Cannot find solution: ", reason);
            }

          El::Matrix<El::BigFloat> y(dual_objective_b_star.Height(), 1);
          El::DistMatrix<El::BigFloat, El::STAR, El::STAR> yp_star(
            solver.y.blocks.at(0));

          El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                   yp_to_y_star.LockedMatrix(), yp_star.LockedMatrix(),
                   El::BigFloat(0.0), y);

          fill_weights(y, max_index, normalization, weights);
          if(rank == 0 && parameters_in.verbosity >= Verbosity::regular)
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
          find_new_points(num_blocks, rank, num_procs,
                          parameters_in.mesh_threshold, epsilon, infinity,
                          function_blocks, weights, points, new_points,
                          has_new_points);
          if(!has_new_points)
            {
              if(parameters.solver.duality_gap_threshold
                 == parameters_in.solver.duality_gap_threshold)
                {
                  parameters.solver.duality_gap_threshold = 0;
                }
              else
                {
                  parameters.solver.duality_gap_threshold
                    = std::max(parameters.solver.duality_gap_threshold
                                 / parameters.duality_gap_reduction,
                               parameters_in.solver.duality_gap_threshold);
                }
            }
        }
      El::DistMatrix<El::BigFloat, El::STAR, El::STAR> yp_star(
        solver.y.blocks.front());
      copy_matrix(yp_star, yp_saved);
      save_checkpoint(parameters.solver.checkpoint_out,
                      parameters_in.verbosity, yp_to_y_star,
                      dual_objective_b_star, yp_saved, points, infinity,
                      parameters.solver.duality_gap_threshold, primal_c_scale,
                      backup_generation, current_generation);
    }
  return weights;
}
