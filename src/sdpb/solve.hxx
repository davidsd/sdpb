#pragma once

#include "SDPB_Parameters.hxx"
#include "save_solution.hxx"
#include "write_timing.hxx"
#include "sdp_solve/memory_estimates.hxx"
#include "sdp_solve/SDP_Solver/run/constraint_matrix_weighted_sum.hxx"
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

// sdpb --reconstructXYmu mu: after a text checkpoint with only x and y,
// X = sum_p A_p x_p and Y = mu X^{-1}. Without the option, such a
// checkpoint is refused: X = Y = 0 is no starting point.
inline void
reconstruct_XY_mu(const SDPB_Parameters &parameters,
                  const Block_Info &block_info, const SDP &sdp,
                  SDP_Solver &solver, Timers &timers)
{
  const auto &ck = parameters.solver.checkpoint_in;
  if(parameters.reconstruct_XY_mu.empty())
    {
      ASSERT(!solver.loaded_without_XY, "Text checkpoint ", ck,
             " has no X_matrix/Y_matrix files; use --reconstructXYmu.");
      return;
    }
  ASSERT(solver.loaded_without_XY,
         "--reconstructXYmu needs a text checkpoint with x_*.txt, y.txt and "
         "no X/Y matrices: ",
         ck);
  Scoped_Timer reconstruct_timer(timers, "reconstruct_XY");
  const El::BigFloat mu(parameters.reconstruct_XY_mu);
  ASSERT(mu > 0, "--reconstructXYmu must be positive: ",
         parameters.reconstruct_XY_mu);
  // X = sum_p A_p x_p  (C = 0, see compute_primal_residues_and_error_P_Ax_X)
  constraint_matrix_weighted_sum(block_info, sdp, solver.x, solver.X);
  // Y = mu X^{-1} = mu L^{-T} L^{-1}, X = L L^T
  for(size_t b = 0; b < solver.X.blocks.size(); ++b)
    {
      auto &X_block = solver.X.blocks[b];
      auto &Y_block = solver.Y.blocks[b];
      if(X_block.Height() == 0)
        continue;
      El::DistMatrix<El::BigFloat> L(X_block);
      try
        {
          El::Cholesky(El::UpperOrLowerNS::LOWER, L);
        }
      catch(std::exception &e)
        {
          RUNTIME_ERROR("--reconstructXYmu: X = sum_p A_p x_p is not "
                        "positive definite in local block ",
                        b, " (", e.what(),
                        "); x must be primal feasible, or nearly.");
        }
      El::DistMatrix<El::BigFloat> L_inv(L.Grid());
      El::Identity(L_inv, L.Height(), L.Width());
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L, L_inv);
      El::Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE, mu,
               L_inv, El::BigFloat(0), Y_block);
      El::MakeSymmetric(El::UpperOrLowerNS::LOWER, Y_block);
    }
  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    std::cout << "Reconstructed X = sum_p A_p x_p and Y = mu X^{-1}, mu = "
              << parameters.reconstruct_XY_mu << '\n';
}

inline void reconstruct_XY_mu(const SDPB_Parameters &parameters,
                              const Sdpb::Sdpa::Block_Info &,
                              const Sdpb::Sdpa::SDP &,
                              Sdpb::Sdpa::SDP_Solver &, Timers &)
{
  ASSERT(parameters.reconstruct_XY_mu.empty(),
         "--reconstructXYmu is implemented for polynomial matrix programs "
         "only, not for SDPA input.");
}

template <class TSolver>
SDP_Solver_Terminate_Reason
solve(const typename TSolver::Block_Info_Type &block_info,
      const SDPB_Parameters &parameters, const Environment &env,
      const std::chrono::time_point<std::chrono::high_resolution_clock>
        &start_time,
      Timers &timers, El::Matrix<int32_t> &block_timings_ms)
{
  Scoped_Timer solve_timer(timers, "sdpb.solve");

  El::Grid grid(block_info.mpi_comm.value);

  auto sdp = read_sdp(parameters, block_info, grid, timers);
  if(parameters.verbosity >= Verbosity::debug)
    {
      print_allocation_message_per_node(env, "SDP", get_allocated_bytes(sdp));
    }

  auto solver = create_solver(parameters, block_info, sdp, grid, timers);
  reconstruct_XY_mu(parameters, block_info, sdp, solver, timers);
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

  if(parameters.verbosity >= Verbosity::debug)
    {
      try
        {
          write_profiling(parameters.solver.checkpoint_out, timers);
        }
      catch(std::exception &e)
        {
          El::Output("An exception has been thrown in write_profiling():");
          El::ReportException(e);
        }
    }

  if(reason == SDP_Solver_Terminate_Reason::SIGTERM_Received)
    {
      if(El::mpi::Rank() == 0)
        El::Output("Received SIGTERM, exiting gracefully...");
      MPI_Finalize();
      std::exit(SIGTERM);
    }
  return reason;
}
