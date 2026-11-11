//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "SDPB_Parameters.hxx"
#include "time_and_solve.hxx"
#include "resident/Resident.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdp_solve/SDP_Solver/run/multimodular/multimodular_blas.hxx"
#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpb_version.hxx"
#include "sdpb_util/Proc_Meminfo.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <El.hpp>

void run_resident(const Environment &env, const SDPB_Parameters &parameters,
                  const std::vector<std::string> &startup_args,
                  mp_bitcnt_t startup_precision, const Block_Info &block_info);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      // The options are parsed at GMP's precision before set_precision()
      // (resident mode parses its solve commands at the same precision)
      const mp_bitcnt_t startup_precision = El::gmp::Precision();
      SDPB_Parameters parameters(argc, argv, env);
      if(!parameters.is_valid())
        {
          return 0;
        }

#ifndef SDPB_GPU
      ASSERT(!parameters.solver.gpu,
             "--gpu: this sdpb was built without the GPU backend "
             "(./waf configure --enable-gpu).");
#endif
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

      if(parameters.is_sdpa())
        {
          ASSERT(!parameters.solver.gpu,
                 "--gpu is implemented for polynomial matrix programs only, "
                 "not for SDPA input.");
          ASSERT(parameters.control_dir.empty(),
                 "--controlDir (resident mode) is implemented for polynomial "
                 "matrix programs only, not for SDPA input.");
          if(parameters.solver.multimodular_trsm)
            RUNTIME_ERROR("--multimodularTrsm is not supported for SDPA "
                          "input (.dat, .dat-s).");
          time_and_solve<Sdpb::Sdpa::SDP_Solver>(env, parameters, start_time);
        }
      else if(!parameters.control_dir.empty())
        {
          if(parameters.solver.multimodular_trsm)
            check_multimodular_precision(Multimodular_Blas_Config(),
                                         El::gmp::Precision());
          // Resident mode (sdpb --controlDir): after the timing run, solve
          // on command
          const std::vector<std::string> startup_args(argv + 1, argv + argc);
          time_and_solve<SDP_Solver>(
            env, parameters, start_time,
            [&](const Block_Info &block_info, const SDPB_Parameters &p) {
              run_resident(env, p, startup_args, startup_precision,
                           block_info);
            });
        }
      else
        {
          if(parameters.solver.multimodular_trsm)
            check_multimodular_precision(Multimodular_Blas_Config(),
                                         El::gmp::Precision());
          time_and_solve<SDP_Solver>(env, parameters, start_time);
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
