//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "SDP_Solver_Parameters.hxx"
#include "Block_Info.hxx"
#include "../Timers.hxx"

#include <El.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

Timers
solve(const boost::filesystem::path &sdp_directory,
      const boost::filesystem::path &out_file,
      const boost::filesystem::path &checkpoint_in,
      const boost::filesystem::path &checkpoint_out,
      const Block_Info &block_info, const SDP_Solver_Parameters &parameters);

void write_timing(const boost::filesystem::path &checkpoint_out,
                  const Block_Info &block_info, const Timers &timers,
                  const bool &debug);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      boost::filesystem::path sdp_directory, out_file, checkpoint_in,
        checkpoint_out, param_file;
      int verbosity;
      SDP_Solver_Parameters parameters;

      po::options_description basic_options("Basic options");
      basic_options.add_options()("help,h", "Show this helpful message.")(
        "sdpDir,s",
        po::value<boost::filesystem::path>(&sdp_directory)->required(),
        "Directory containing preprocessed SDP data files.")(
        "paramFile,p", po::value<boost::filesystem::path>(&param_file),
        "Any parameter can optionally be set via this file in key=value "
        "format. Command line arguments override values in the parameter "
        "file.")("outFile,o", po::value<boost::filesystem::path>(&out_file),
                 "The optimal solution is saved to this file in Mathematica "
                 "format. Defaults to sdpDir with '.out' extension.")(
        "checkpointDir,c", po::value<boost::filesystem::path>(&checkpoint_out),
        "Checkpoints are saved to this directory every checkpointInterval. "
        "Defaults to sdpDir with '.ck' extension.")(
        "initialCheckpointDir,i",
        po::value<boost::filesystem::path>(&checkpoint_in),
        "The initial checkpoint directory to load. Defaults to "
        "checkpointDir.")("verbosity",
                          po::value<int>(&verbosity)->default_value(1),
                          "Verbosity.  0 -> no output, 1 -> regular "
                          "output, 2 -> debug output");

      // We set default parameters using El::BigFloat("1e-10",10)
      // rather than a straight double precision 1e-10 so that results
      // are more reproducible at high precision.  Using double
      // precision defaults results in differences of about 1e-15 in
      // primalObjective after one step.
      po::options_description solver_params_options("Solver parameters");
      solver_params_options.add_options()(
        "precision", po::value<int>(&parameters.precision)->default_value(400),
        "The precision, in the number of bits, for numbers in the "
        "computation. "
        " This should be less than or equal to the precision used when "
        "preprocessing the XML input files with 'pvm2sdp'.  GMP will round "
        "this up to a multiple of 32 or 64, depending on the system.")(
        "checkpointInterval",
        po::value<int>(&parameters.checkpoint_interval)->default_value(3600),
        "Save checkpoints to checkpointDir every checkpointInterval "
        "seconds.")(
        "noFinalCheckpoint",
        po::bool_switch(&parameters.no_final_checkpoint)->default_value(false),
        "Don't save a final checkpoint after terminating (useful when "
        "debugging).")(
        "findPrimalFeasible",
        po::bool_switch(&parameters.find_primal_feasible)->default_value(false),
        "Terminate once a primal feasible solution is found.")(
        "findDualFeasible",
        po::bool_switch(&parameters.find_dual_feasible)->default_value(false),
        "Terminate once a dual feasible solution is found.")(
        "detectPrimalFeasibleJump",
        po::bool_switch(&parameters.detect_primal_feasible_jump)
          ->default_value(false),
        "Terminate if a primal-step of 1 is taken. This often indicates that "
        "a "
        "primal feasible solution would be found if the precision were high "
        "enough. Try increasing either primalErrorThreshold or precision "
        "and run from the latest checkpoint.")(
        "detectDualFeasibleJump",
        po::bool_switch(&parameters.detect_dual_feasible_jump)
          ->default_value(false),
        "Terminate if a dual-step of 1 is taken. This often indicates that a "
        "dual feasible solution would be found if the precision were high "
        "enough. Try increasing either dualErrorThreshold or precision "
        "and run from the latest checkpoint.")(
        "maxIterations",
        po::value<int>(&parameters.max_iterations)->default_value(500),
        "Maximum number of iterations to run the solver.")(
        "maxRuntime",
        po::value<int>(&parameters.max_runtime)->default_value(86400),
        "Maximum amount of time to run the solver in seconds.")(
        "procsPerNode", po::value<int>(&parameters.procs_per_node)->required(),
        "This option is **required**.\n\n"
        "The number of processes that can run on a node.  When running on "
        "more "
        "than one node, the load balancer needs to know how many processes "
        "are assigned to each node.  On a laptop or desktop, this would be "
        "the number of physical cores on your machine, not including "
        "hyperthreaded cores.  For current laptops (2018), this is probably "
        "2 or 4.\n\n"
        "If you are using the Slurm workload manager, this should be set to "
        "'$SLURM_NTASKS_PER_NODE'.")(
        "dualityGapThreshold",
        po::value<El::BigFloat>(&parameters.duality_gap_threshold)
          ->default_value(El::BigFloat("1e-30", 10)),
        "Threshold for duality gap (roughly the difference in primal and dual "
        "objective) at which the solution is considered "
        "optimal. Corresponds to SDPA's epsilonStar.")(
        "primalErrorThreshold",
        po::value<El::BigFloat>(&parameters.primal_error_threshold)
          ->default_value(El::BigFloat("1e-30", 10)),
        "Threshold for feasibility of the primal problem. Corresponds to "
        "SDPA's epsilonBar.")(
        "dualErrorThreshold",
        po::value<El::BigFloat>(&parameters.dual_error_threshold)
          ->default_value(El::BigFloat("1e-30", 10)),
        "Threshold for feasibility of the dual problem. Corresponds to SDPA's "
        "epsilonBar.")(
        "initialMatrixScalePrimal",
        po::value<El::BigFloat>(&parameters.initial_matrix_scale_primal)
          ->default_value(El::BigFloat("1e20", 10)),
        "The primal matrix X begins at initialMatrixScalePrimal times the "
        "identity matrix. Corresponds to SDPA's lambdaStar.")(
        "initialMatrixScaleDual",
        po::value<El::BigFloat>(&parameters.initial_matrix_scale_dual)
          ->default_value(El::BigFloat("1e20", 10)),
        "The dual matrix Y begins at initialMatrixScaleDual times the "
        "identity matrix. Corresponds to SDPA's lambdaStar.")(
        "feasibleCenteringParameter",
        po::value<El::BigFloat>(&parameters.feasible_centering_parameter)
          ->default_value(El::BigFloat("0.1", 10)),
        "Shrink the complementarity X Y by this factor when the primal and "
        "dual "
        "problems are feasible. Corresponds to SDPA's betaStar.")(
        "infeasibleCenteringParameter",
        po::value<El::BigFloat>(&parameters.infeasible_centering_parameter)
          ->default_value(El::BigFloat("0.3", 10)),
        "Shrink the complementarity X Y by this factor when either the primal "
        "or dual problems are infeasible. Corresponds to SDPA's betaBar.")(
        "stepLengthReduction",
        po::value<El::BigFloat>(&parameters.step_length_reduction)
          ->default_value(El::BigFloat("0.7", 10)),
        "Shrink each newton step by this factor (smaller means slower, more "
        "stable convergence). Corresponds to SDPA's gammaStar.")(
        "maxComplementarity",
        po::value<El::BigFloat>(&parameters.max_complementarity)
          ->default_value(El::BigFloat("1e100", 10)),
        "Terminate if the complementarity mu = Tr(X Y)/dim(X) "
        "exceeds this value.");

      po::options_description cmd_line_options;
      cmd_line_options.add(basic_options).add(solver_params_options);

      po::variables_map variables_map;

      try
        {
          po::store(po::parse_command_line(argc, argv, cmd_line_options),
                    variables_map);

          if(variables_map.count("help") != 0)
            {
              if(El::mpi::Rank() == 0)
                {
                  std::cout << cmd_line_options << '\n';
                }
              return 0;
            }

          if(variables_map.count("paramFile") != 0)
            {
              param_file
                = variables_map["paramFile"].as<boost::filesystem::path>();
              std::ifstream ifs(param_file.string().c_str());
              if(!ifs.good())
                {
                  throw std::runtime_error("Could not open '"
                                           + param_file.string() + "'");
                }

              po::store(po::parse_config_file(ifs, solver_params_options),
                        variables_map);
            }

          po::notify(variables_map);

          if(!boost::filesystem::exists(sdp_directory))
            {
              throw std::runtime_error("sdp directory '"
                                       + sdp_directory.string()
                                       + "' does not exist");
            }
          if(!boost::filesystem::is_directory(sdp_directory))
            {
              throw std::runtime_error("sdp directory '"
                                       + sdp_directory.string()
                                       + "' is not a directory");
            }

          if(variables_map.count("outFile") == 0)
            {
              out_file = sdp_directory;
              if(out_file.filename() == ".")
                {
                  out_file = out_file.parent_path();
                }
              out_file.replace_extension("out");
            }

          if(variables_map.count("checkpointDir") == 0)
            {
              checkpoint_out = sdp_directory;
              if(checkpoint_out.filename() == ".")
                {
                  checkpoint_out = checkpoint_out.parent_path();
                }
              checkpoint_out.replace_extension("ck");
            }

          if(variables_map.count("initialCheckpointDir") == 0)
            {
              checkpoint_in = checkpoint_out;
            }

          if(El::mpi::Rank() == 0)
            {
              std::ofstream ofs(out_file.string().c_str());
              if(!ofs)
                {
                  std::cerr << "Cannot write to outFile." << '\n';
                  return 1;
                }
            }

          if(verbosity != 0 && verbosity != 1 && verbosity != 2)
            {
              std::cerr << "Invalid number for Verbosity.  Only 0, 1 or 2 are "
                           "allowed\n";
              return 1;
            }
          else
            {
              parameters.verbosity = static_cast<Verbosity>(verbosity);
            }
        }
      catch(po::error &e)
        {
          El::ReportException(e);
          El::mpi::Abort(El::mpi::COMM_WORLD, 1);
        }

      El::gmp::SetPrecision(parameters.precision);
      if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          std::cout << "SDPB started at "
                    << boost::posix_time::second_clock::local_time() << '\n'
                    << "SDP directory   : " << sdp_directory << '\n'
                    << "out file        : " << out_file << '\n'
                    << "checkpoint in   : " << checkpoint_in << '\n'
                    << "checkpoint out  : " << checkpoint_out << '\n'
                    << "\nParameters:\n"
                    << parameters << '\n';
        }

      Block_Info block_info(sdp_directory, checkpoint_in,
                            parameters.procs_per_node, parameters.verbosity);
      // Only generate a block_timings file if
      // 1) We are running in parallel
      // 2) We did not load a block_timings file
      // 3) We are not going to load a checkpoint.
      if(El::mpi::Size(El::mpi::COMM_WORLD) > 1
         && block_info.block_timings_filename.empty()
         && !exists(checkpoint_in
                    / ("checkpoint." + std::to_string(El::mpi::Rank()))))
        {
          if(parameters.verbosity >= Verbosity::regular
             && El::mpi::Rank() == 0)
            {
              std::cout << "Performing a timing run\n";
            }
          SDP_Solver_Parameters timing_parameters(parameters);
          timing_parameters.max_iterations = 2;
          timing_parameters.no_final_checkpoint = true;
          timing_parameters.checkpoint_interval=std::numeric_limits<int>::max();
          timing_parameters.verbosity = Verbosity::none;
          Timers timers(solve(sdp_directory, out_file, checkpoint_in,
                              checkpoint_out, block_info, timing_parameters));

          write_timing(checkpoint_out, block_info, timers,
                       timing_parameters.verbosity >= Verbosity::debug);
          El::mpi::Barrier(El::mpi::COMM_WORLD);
          block_info
            = Block_Info(sdp_directory, checkpoint_out,
                         parameters.procs_per_node, parameters.verbosity);

          parameters.max_runtime-=timers.front().second.elapsed_seconds();
        }
      else if(!block_info.block_timings_filename.empty()
              && block_info.block_timings_filename
                   != (checkpoint_out / "block_timings"))
        {
          if(El::mpi::Rank() == 0)
            {
              create_directories(checkpoint_out);
              copy_file(block_info.block_timings_filename,
                        checkpoint_out / "block_timings");
            }
        }
      solve(sdp_directory, out_file, checkpoint_in, checkpoint_out, block_info,
            parameters);
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
