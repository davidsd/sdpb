#pragma once

#include "create_block_info.hxx"
#include "solve.hxx"
#include "write_timing.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdpb_util/Environment.hxx"
#include "SDPB_Parameters.hxx"
#include "sdpb_util/malloc_trim.hxx"

#include <boost/date_time/posix_time/posix_time.hpp>

#include <chrono>
#include <filesystem>

template <class TSolver>
Timers time_and_solve(
  const Environment &env, SDPB_Parameters parameters,
  const std::chrono::time_point<std::chrono::high_resolution_clock> &start_time)
{
  namespace fs = std::filesystem;
  const auto block_timings_path = get_block_timings_path(parameters);
  auto block_info
    = create_block_info<TSolver>(env, parameters, block_timings_path);
  // Only generate a block_timings file if
  // 1) We are running in parallel
  // 2) We did not load a block_timings file
  // 3) We are not going to load a checkpoint.
  if(El::mpi::Size(El::mpi::COMM_WORLD) > 1 && block_timings_path.empty()
     && !exists(parameters.solver.checkpoint_in / "checkpoint.0"))
    {
      if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          El::Output(boost::posix_time::second_clock::local_time(),
                     " Start timing run");
        }
      SDPB_Parameters timing_parameters(parameters);
      timing_parameters.solver.max_iterations = 2;
      timing_parameters.no_final_checkpoint = true;
      timing_parameters.solver.checkpoint_interval
        = std::numeric_limits<int64_t>::max();
      timing_parameters.solver.max_runtime
        = std::numeric_limits<int64_t>::max();
      timing_parameters.solver.duality_gap_threshold = 0;
      timing_parameters.solver.primal_error_threshold = 0;
      timing_parameters.solver.dual_error_threshold = 0;
      timing_parameters.solver.min_primal_step = 0;
      timing_parameters.solver.min_dual_step = 0;
      if(timing_parameters.verbosity < Verbosity::debug)
        {
          timing_parameters.verbosity = Verbosity::none;
        }
      El::Matrix<int32_t> block_timings_ms;
      Timers timers = solve<TSolver>(block_info, timing_parameters, env,
                                     start_time, block_timings_ms);
      // Release memory back to OS (glibc-only)
      malloc_trim(env, timers);

      if(block_timings_ms.Height() == 0 && block_timings_ms.Width() == 0)
        {
          RUNTIME_ERROR(
            "block_timings vector is empty, probably because timing run "
            "exited before completing two solver iterations.");
        }

      write_block_timings(timing_parameters.solver.checkpoint_out,
                          block_timings_ms, timing_parameters.verbosity);
      if(timing_parameters.verbosity >= Verbosity::debug)
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

      El::mpi::Barrier(El::mpi::COMM_WORLD);
      auto new_info
        = create_block_info<TSolver>(env, parameters, block_timings_ms);
      swap(block_info, new_info);

      const auto elapsed_seconds
        = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::high_resolution_clock::now() - start_time)
            .count();
      parameters.solver.max_runtime -= elapsed_seconds;
    }
  else if(!block_timings_path.empty()
          && block_timings_path
               != (parameters.solver.checkpoint_out / "block_timings"))
    {
      if(El::mpi::Rank() == 0)
        {
          create_directories(parameters.solver.checkpoint_out);
          copy_file(block_timings_path,
                    parameters.solver.checkpoint_out / "block_timings",
                    fs::copy_options::overwrite_existing);
        }
    }
  El::Matrix<int32_t> block_timings_ms;
  return solve<TSolver>(block_info, parameters, env, start_time,
                        block_timings_ms);
}
