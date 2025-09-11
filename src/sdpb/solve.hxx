#pragma once

#include "SDPB_Parameters.hxx"
#include "save_solution.hxx"
#include "sdp_solve/memory_estimates.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdpa_solve/memory_estimates.hxx"
#include "sdpa_solve/sdpa_solve.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>
#include <csignal>
#include <filesystem>
#include <boost/date_time/posix_time/posix_time.hpp>

inline SDP
read_sdp(const SDPB_Parameters &parameters, const Block_Info &block_info,
         const El::Grid &grid, Timers &timers)
{
  Scoped_Timer read_sdp_timer(timers, "read_sdp");
  SDP sdp(parameters.sdp_path, block_info, grid, timers);
  if(El::mpi::Rank() == 0 && parameters.write_solution.vector_z)
    {
      ASSERT(sdp.normalization.has_value(),
             "Please provide SDP with valid normalization.json "
             "or exclude z from --writeSolution arguments.");
    }
  return sdp;
}

Sdpb::Sdpa::SDP read_sdp(const SDPB_Parameters &parameters,
                         const Sdpb::Sdpa::Block_Info &block_info,
                         const El::Grid &grid, Timers &timers)
{
  Scoped_Timer read_sdp_timer(timers, "read_sdp");
  return Sdpb::Sdpa::SDP(parameters.sdp_path, block_info, grid, timers);
}

inline SDP_Solver
create_solver(const SDPB_Parameters &parameters, const Block_Info &block_info,
              const SDP &sdp, const El::Grid &grid, Timers &timers)
{
  Scoped_Timer solver_ctor_timer(timers, "SDP_Solver.ctor");
  return SDP_Solver(parameters.solver, parameters.verbosity,
                    parameters.require_initial_checkpoint, block_info, grid,
                    sdp.dual_objective_b.Height());
}

inline Sdpb::Sdpa::SDP_Solver
create_solver(const SDPB_Parameters &parameters,
              const Sdpb::Sdpa::Block_Info &block_info,
              [[maybe_unused]] const Sdpb::Sdpa::SDP &, const El::Grid &grid,
              Timers &timers)
{
  Scoped_Timer solver_ctor_timer(timers, "SDP_Solver.ctor");
  return Sdpb::Sdpa::SDP_Solver(parameters.solver, parameters.verbosity,
                                parameters.require_initial_checkpoint,
                                block_info, grid);
}

template <class TSolver>
Timers solve(const typename TSolver::Block_Info_Type &block_info,
             const SDPB_Parameters &parameters, const Environment &env,
             const std::chrono::time_point<std::chrono::high_resolution_clock>
               &start_time,
             El::Matrix<int32_t> &block_timings_ms)
{
  Timers timers(env, parameters.verbosity);
  Scoped_Timer solve_timer(timers, "sdpb.solve");

  El::Grid grid(block_info.mpi_comm.value);

  auto sdp = read_sdp(parameters, block_info, grid, timers);
  if(parameters.verbosity >= Verbosity::debug)
    {
      print_allocation_message_per_node(env, "SDP", get_allocated_bytes(sdp));
    }

  auto solver = create_solver(parameters, block_info, sdp, grid, timers);
  if(parameters.verbosity >= Verbosity::debug)
    {
      print_allocation_message_per_node(env, "SDP_Solver",
                                        get_allocated_bytes(solver));
    }

  const boost::property_tree::ptree parameters_tree(
    to_property_tree(parameters));

  const auto iterations_json_path
    = parameters.out_directory / "iterations.json";
  SDP_Solver_Terminate_Reason reason(solver.run(
    env, parameters.solver, parameters.verbosity, parameters_tree, block_info,
    sdp, grid, start_time, iterations_json_path, timers, block_timings_ms));

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      set_stream_precision(std::cout);
      std::cout << "-----" << reason << "-----\n"
                << '\n'
                << "primalObjective = " << solver.primal_objective << '\n'
                << "dualObjective   = " << solver.dual_objective << '\n'
                << "dualityGap      = " << solver.duality_gap << '\n'
                << "primalError     = " << primal_error(solver) << '\n'
                << "dualError       = " << solver.dual_error << '\n'
                << '\n';
    }

  if(reason == SDP_Solver_Terminate_Reason::SIGTERM_Received
     || !parameters.no_final_checkpoint)
    {
      Scoped_Timer save_timer(timers, "save_checkpoint");
      solver.save_checkpoint(parameters.solver.checkpoint_out,
                             parameters.verbosity, parameters_tree);
    }

  {
    Scoped_Timer save_timer(timers, "save_solution");
    const auto runtime = std::chrono::duration_cast<std::chrono::seconds>(
                     std::chrono::high_resolution_clock::now() - start_time)
                     .count();
    save_solution(solver, block_info, sdp, parameters, reason, runtime);
  }

  if(reason == SDP_Solver_Terminate_Reason::SIGTERM_Received)
    {
      if(El::mpi::Rank() == 0)
        El::Output("Received SIGTERM, exiting gracefully...");
      MPI_Finalize();
      std::exit(SIGTERM);
    }
  return timers;
}
