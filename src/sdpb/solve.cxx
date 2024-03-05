#include "SDPB_Parameters.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <filesystem>

namespace fs = std::filesystem;

void save_solution(
  const SDP_Solver &solver,
  const SDP_Solver_Terminate_Reason &terminate_reason, const int64_t &runtime,
  const fs::path &out_directory, const Write_Solution &write_solution,
  const std::vector<size_t> &block_indices,
  const std::optional<std::vector<El::BigFloat>> &normalization,
  const Verbosity &verbosity);

Timers solve(const Block_Info &block_info, const SDPB_Parameters &parameters,
             const Environment &env,
             const std::chrono::time_point<std::chrono::high_resolution_clock>
               &start_time)
{
  Timers timers(env, parameters.verbosity >= Verbosity::debug);
  Scoped_Timer solve_timer(timers, "sdpb.solve");

  El::Grid grid(block_info.mpi_comm.value);

  Scoped_Timer read_sdp_timer(timers, "read_sdp");
  SDP sdp(parameters.sdp_path, block_info, grid, timers);
  if(El::mpi::Rank() == 0 && parameters.write_solution.vector_z)
    {
      ASSERT(sdp.normalization.has_value(),
             "Please provide SDP with valid normalization.json "
             "or exclude z from --writeSolution arguments.");
    }
  read_sdp_timer.stop();

  Scoped_Timer solver_ctor_timer(timers, "SDP_Solver.ctor");
  SDP_Solver solver(parameters.solver, parameters.verbosity,
                    parameters.require_initial_checkpoint, block_info, grid,
                    sdp.dual_objective_b.Height());
  solver_ctor_timer.stop();

  const boost::property_tree::ptree parameters_tree(
    to_property_tree(parameters));
  SDP_Solver_Terminate_Reason reason(
    solver.run(parameters.solver, parameters.verbosity, parameters_tree,
               block_info, sdp, grid, start_time, timers));

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
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

  if(!parameters.no_final_checkpoint)
    {
      solver.save_checkpoint(parameters.solver.checkpoint_out,
                             parameters.verbosity, parameters_tree);
    }

  auto runtime = std::chrono::duration_cast<std::chrono::seconds>(
                   std::chrono::high_resolution_clock::now() - start_time)
                   .count();
  save_solution(solver, reason, runtime, parameters.out_directory,
                parameters.write_solution, block_info.block_indices,
                sdp.normalization, parameters.verbosity);
  return timers;
}
