//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "SDPB_Parameters.hxx"
#include "time_and_solve.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpb_util/Proc_Meminfo.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <El.hpp>

#include <filesystem>

namespace fs = std::filesystem;

bool is_sdpa_file(fs::path sdp_path)
{
  if(sdp_path.extension() == ".gz")
    sdp_path.replace_extension();
  if(sdp_path.extension() == ".dat" || sdp_path.extension() == ".dat-s")
    return true;
  return false;
}

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      SDPB_Parameters parameters(argc, argv, env);
      if(!parameters.is_valid())
        {
          return 0;
        }

      Environment::set_precision(parameters.solver.precision);
      auto start_time = std::chrono::high_resolution_clock::now();
      if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          // Print command line
          if(parameters.verbosity >= Verbosity::debug)
            {
              std::vector<std::string> arg_list(argv, argv + argc);
              for(const auto &arg : arg_list)
                std::cout << arg << " ";
              std::cout << std::endl;
            }
          std::cout << boost::posix_time::second_clock::local_time()
                    << " Start SDPB" << '\n'
                    << "SDPB version: " << SDPB_VERSION_STRING << '\n'
                    << "MPI processes: " << El::mpi::Size()
                    << ", nodes: " << env.num_nodes() << '\n'
                    << parameters << std::endl;
        }

      if(parameters.verbosity >= Verbosity::debug)
        {
          if(env.comm_shared_mem.Rank() == 0)
            {
              bool res;
              auto meminfo = Proc_Meminfo::try_read(res, true);
              if(res)
                {
                  El::Output("node=", env.node_index(), ": MemUsed: ",
                             pretty_print_bytes(meminfo.mem_used()));
                }
            }
          // Make sure that we don't allocate anything before printing MemUsed
          El::mpi::Barrier(env.comm_shared_mem);
        }

      const Timers timers
        = is_sdpa_file(parameters.sdp_path)
            ? time_and_solve<Sdpb::Sdpa::SDP_Solver>(env, parameters,
                                                     start_time)
            : time_and_solve<SDP_Solver>(env, parameters, start_time);
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
