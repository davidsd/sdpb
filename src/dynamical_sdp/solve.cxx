#include "Dynamical_Parameters.hxx"
#include "../dynamical_solve.hxx"
#include "../set_stream_precision.hxx"

#include <El.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>

void save_solution(const Dynamical_Solver &solver, const Dynamical_Solver_Terminate_Reason,
                   const std::pair<std::string, Timer> &timer_pair,
                   const boost::filesystem::path &out_directory,
                   const Write_Solution &write_solution,
                   const std::vector<size_t> &block_indices,
                   const Verbosity &verbosity, const El::Matrix<El::BigFloat> &extParamStep);

Timers solve(const Block_Info &block_info, const Dynamical_Parameters &parameters)
{
  El::Grid grid(block_info.mpi_comm.value);
  SDP sdp(parameters.sdp_path, block_info, grid);
  Dynamical_Solver solver(parameters.solver, parameters.verbosity,
                    parameters.require_initial_checkpoint, block_info, grid,
                    sdp.dual_objective_b.Height());
  bool update_sdp = false; 

  Timers timers(parameters.verbosity >= Verbosity::debug);
  const boost::property_tree::ptree parameters_tree(
    to_property_tree(parameters));
  El::Matrix<El::BigFloat> extParamStep;

  const boost::filesystem::path out_path(parameters.out_directory / "externalParamStep.txt");
  Dynamical_Solver_Terminate_Reason reason(solver.run_dynamical(
    parameters.solver, parameters.verbosity,
    sdp, parameters_tree, block_info, grid, timers, update_sdp, extParamStep));

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      set_stream_precision(std::cout);
	  std::cout << "-----" << reason << "-----\n"
		  << '\n'
		  << "primalObjective        = " << solver.primal_objective << '\n'
		  << "dualObjective          = " << solver.dual_objective << '\n'
		  << "dualityGap             = " << solver.duality_gap << '\n'
		  << "primalError            = " << solver.primal_error() << '\n'
		  << "dualError              = " << solver.dual_error << '\n'
		  << "extStepSize            = " << solver.external_step_size << '\n'
		  << "totalIteration so far  = " << solver.total_iteration << '\n'
		  << "BFGSHessianUpdated     = " << solver.hess_BFGS_updateQ << '\n'
		  << '\n';
    }

  if(!parameters.no_final_checkpoint)
    {
      solver.save_checkpoint(parameters.solver.solver_parameters.checkpoint_out, parameters.verbosity,
                             parameters_tree);
    }
  save_solution(solver, reason, timers.front(), parameters.out_directory,
                parameters.write_solution, block_info.block_indices,
                parameters.verbosity, extParamStep);
  return timers;
}
