//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "Dynamical_Parameters.hxx"

#include <El.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace fs = std::filesystem;

Timers solve(const Block_Info &block_info,
             const Dynamical_Parameters &parameters, const Environment &env,
             const std::chrono::time_point<std::chrono::high_resolution_clock>
               &start_time,
             El::Matrix<int32_t> &block_timings_ms);

void write_block_timings(const fs::path &checkpoint_out,
                         const Block_Info &block_info,
                         const El::Matrix<int32_t> &block_timings_ms,
                         const Verbosity verbosity);

void write_profiling(const fs::path &checkpoint_out, const Timers &timers);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  PrecParameters prec_parameter(argc, argv);

  if(El::mpi::Rank() == 0)
    std::cout << "Elemental precision = " << prec_parameter.prec << "\n"
              << std::flush;

  Environment::set_precision(prec_parameter.prec);

  /*
  // a very strange bug : if I compile with boost 1.77 , if I also use "mpirun -n 2", then somehow the copy stuck 
  // at copy_file, even if I re-run it with "mpirun -n 1"
  // compile with boost 1.67 is fine

  std::cout << "before copy_file\n" << std::flush;

  if (El::mpi::Rank() == 0)
  {
	  copy_file(
		  "./Proj_Ising_pd11_3d_pureBFGS_newstrategy_debug_mpi_ck/SDPBFiles/0.519000000000_1.40000000000_0.576792155454_0.816890940950_Aug30_08h16m15s.ck/"
		  "checkpoint.json",
		  "./Proj_Ising_pd11_3d_pureBFGS_newstrategy_debug_mpi_ck/SDPBFiles/0.519218317386_1.40203144757_0.575973071057_0.817468666933_Aug30_08h48m54s.ck/"
		  "checkpoint.json",
		  std::filesystem::copy_options::overwrite_existing);
  }

  std::cout << "after copy_file\n" << std::flush;

  return 1;

  */

  try
    {
      Dynamical_Parameters parameters(argc, argv);
      if(!parameters.is_valid())
        {
          return 0;
        }

      //El::gmp::SetPrecision(parameters.solver.solver_parameters.precision);
      auto start_time = std::chrono::high_resolution_clock::now();
      if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          std::cout << boost::posix_time::second_clock::local_time()
                    << " Start Dynamical SDP" << '\n'
                    << "SDPB version: " << SDPB_VERSION_STRING << '\n'
                    << "MPI processes: " << El::mpi::Size()
                    << ", nodes: " << env.num_nodes() << '\n'
                    << parameters << std::endl;
        }

      Block_Info block_info(env, parameters.sdp_path,
                            parameters.solver.solver_parameters.checkpoint_in,
                            parameters.proc_granularity, parameters.verbosity);

      // Only generate a block_timings file if
      // 1) We are running in parallel
      // 2) We did not load a block_timings file
      // 3) We are not going to load a checkpoint.
      if(El::mpi::Size(El::mpi::COMM_WORLD) > 1
         && block_info.block_timings_filename.empty()
         && !exists(parameters.solver.solver_parameters.checkpoint_in
                    / "checkpoint.0"))
        {
          if(parameters.verbosity >= Verbosity::regular
             && El::mpi::Rank() == 0)
            {
              El::Output(boost::posix_time::second_clock::local_time(),
                         " Start timing run");
            }
          Dynamical_Parameters timing_parameters(parameters);
          timing_parameters.solver.solver_parameters.max_iterations = 2;
          timing_parameters.no_final_checkpoint = true;
          timing_parameters.solver.solver_parameters.checkpoint_interval
            = std::numeric_limits<int64_t>::max();
          timing_parameters.solver.solver_parameters.max_runtime
            = std::numeric_limits<int64_t>::max();
          timing_parameters.solver.solver_parameters.duality_gap_threshold = 0;
          timing_parameters.solver.solver_parameters.primal_error_threshold
            = 0;
          timing_parameters.solver.solver_parameters.dual_error_threshold = 0;
          timing_parameters.solver.solver_parameters.min_primal_step = 0;
          timing_parameters.solver.solver_parameters.min_dual_step = 0;
          if(timing_parameters.verbosity != Verbosity::debug)
            {
              timing_parameters.verbosity = Verbosity::none;
            }
          El::Matrix<int32_t> block_timings_ms;
          Timers timers(solve(block_info, timing_parameters, env, start_time,
                              block_timings_ms));

          if(block_timings_ms.Height() == 0 && block_timings_ms.Width() == 0)
            {
              RUNTIME_ERROR(
                "block_timings vector is empty, probably because "
                "timing run exited before completing two solver iterations.");
            }

          write_block_timings(
            timing_parameters.solver.solver_parameters.checkpoint_out,
            block_info, block_timings_ms,
            timing_parameters.verbosity);
          if(timing_parameters.verbosity >= Verbosity::debug)
            {
              try
                {
                  write_profiling(
                    parameters.solver.solver_parameters.checkpoint_out,
                    timers);
                }
              catch(std::exception &e)
                {
                  El::Output(
                    "An exception has been thrown in write_profiling():");
                  El::ReportException(e);
                }
            }

          El::mpi::Barrier(El::mpi::COMM_WORLD);
          Block_Info new_info(env, parameters.sdp_path, block_timings_ms,
                              parameters.proc_granularity,
                              parameters.verbosity);
          std::swap(block_info, new_info);

          auto elapsed_seconds
            = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::high_resolution_clock::now() - start_time)
                .count();
          parameters.solver.solver_parameters.max_runtime -= elapsed_seconds;
        }
      else if(!block_info.block_timings_filename.empty()
              && block_info.block_timings_filename
                   != (parameters.solver.solver_parameters.checkpoint_out
                       / "block_timings"))
        {
          //std::cout << "before copy_file : from " << block_info.block_timings_filename
          //	  << " to " << parameters.solver.solver_parameters.checkpoint_out / "block_timings" << " \n" << std::flush;

          if(El::mpi::Rank() == 0)
            {
              create_directories(
                parameters.solver.solver_parameters.checkpoint_out);

              copy_file(block_info.block_timings_filename,
                        parameters.solver.solver_parameters.checkpoint_out
                          / "block_timings",
                        std::filesystem::copy_options::overwrite_existing);
            }

          //std::cout << "after copy_file \n " << std::flush;
        }
      El::Matrix<int32_t> block_timings_ms;
      Timers timers(
        solve(block_info, parameters, env, start_time, block_timings_ms));
      if(parameters.verbosity >= Verbosity::debug)
        {
          try
            {
              write_profiling(
                parameters.solver.solver_parameters.checkpoint_out, timers);
            }
          catch(std::exception &e)
            {
              El::Output("An exception has been thrown in write_profiling():");
              El::ReportException(e);
            }
        }
    }
  catch(std::exception &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
  catch(...)
    {
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
