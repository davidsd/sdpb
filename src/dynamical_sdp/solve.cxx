#include "Dynamical_Parameters.hxx"
#include "dynamical_solve/dynamical_solve.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>
#include <filesystem>

void save_solution(
  const Dynamical_Solver &solver,
  const Dynamical_Solver_Terminate_Reason &terminate_reason,
  const int64_t &solver_runtime, const std::filesystem::path &out_directory,
  const Write_Solution &write_solution,
  const std::vector<size_t> &block_indices,
  const std::optional<std::vector<El::BigFloat>> &normalization,
  const Verbosity &verbosity, const El::Matrix<El::BigFloat> &extParamStep);

Timers solve(const Block_Info &block_info,
             const Dynamical_Parameters &parameters, const Environment &env,
             const std::chrono::time_point<std::chrono::high_resolution_clock>
               &start_time)
{
  Timers timers(env, parameters.verbosity >= Verbosity::debug);
  Scoped_Timer solve_timer(timers, "dynamical_sdp.solve");

  El::Grid grid(block_info.mpi_comm.value);

  Scoped_Timer read_sdp_timer(timers, "read_sdp");
  SDP sdp(parameters.sdp_path, block_info, grid, timers);
  read_sdp_timer.stop();

  Scoped_Timer solver_ctor_timer(timers, "Dynamical_Solver.ctor");
  Dynamical_Solver solver(parameters.solver, parameters.verbosity,
                          parameters.require_initial_checkpoint, block_info,
                          grid, sdp.dual_objective_b.Height());
  solver_ctor_timer.stop();
  bool update_sdp = false;

  const boost::property_tree::ptree parameters_tree(
    to_property_tree(parameters));
  El::Matrix<El::BigFloat> extParamStep;

  const std::filesystem::path out_path(parameters.out_directory
                                       / "externalParamStep.txt");
  Dynamical_Solver_Terminate_Reason reason(solver.run_dynamical(
    parameters.solver, parameters.verbosity, sdp, parameters_tree, block_info,
    grid, timers, update_sdp, extParamStep));

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      set_stream_precision(std::cout);
      std::cout << "-----" << reason << "-----\n"
                << '\n'
                << "primalObjective        = " << solver.primal_objective
                << '\n'
                << "dualObjective          = " << solver.dual_objective << '\n'
                << "dualityGap             = " << solver.duality_gap << '\n'
                << "primalError            = " << solver.primal_error() << '\n'
                << "dualError              = " << solver.dual_error << '\n'
                << "extStepSize            = " << solver.external_step_size
                << '\n'
                << "totalIteration so far  = " << solver.total_iteration
                << '\n'
                << "BFGSHessianUpdated     = " << solver.hess_BFGS_updateQ
                << '\n'
                << '\n';
    }

  if(!parameters.no_final_checkpoint)
    {
      solver.save_checkpoint(
        parameters.solver.solver_parameters.checkpoint_out,
        parameters.verbosity, parameters_tree);
    }

  auto runtime = std::chrono::duration_cast<std::chrono::seconds>(
                   std::chrono::high_resolution_clock::now() - start_time)
                   .count();
  save_solution(solver, reason, runtime, parameters.out_directory,
                parameters.write_solution, block_info.block_indices,
                sdp.normalization, parameters.verbosity, extParamStep);
  return timers;
}
