#include "../SDP_Solver_Parameters.hxx"

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>

namespace po = boost::program_options;

SDP_Solver_Parameters::SDP_Solver_Parameters(int argc, char *argv[])
{
  int int_verbosity;
  std::string write_solution_string;
  using namespace std::string_literals;

  po::options_description required_options("Required options");
  required_options.add_options()(
    "sdpDir,s", po::value<boost::filesystem::path>(&sdp_directory)->required(),
    "Directory containing preprocessed SDP data files.");
  required_options.add_options()(
    "procsPerNode", po::value<size_t>(&procs_per_node)->required(),
    "The number of processes that can run on a node.  When running on "
    "more "
    "than one node, the load balancer needs to know how many processes "
    "are assigned to each node.  On a laptop or desktop, this would be "
    "the number of physical cores on your machine, not including "
    "hyperthreaded cores.  For current laptops (2018), this is probably "
    "2 or 4.\n\n"
    "If you are using the Slurm workload manager, this should be set to "
    "'$SLURM_NTASKS_PER_NODE'.");

  po::options_description basic_options("Basic options");
  basic_options.add_options()("help,h", "Show this helpful message.");
  basic_options.add_options()("version",
                              "Show version and configuration info.");
  basic_options.add_options()(
    "paramFile,p", po::value<boost::filesystem::path>(&param_file),
    "Any parameter can optionally be set via this file in key=value "
    "format. Command line arguments override values in the parameter "
    "file.");
  basic_options.add_options()(
    "outDir,o", po::value<boost::filesystem::path>(&out_directory),
    "The optimal solution is saved to this directory in Mathematica "
    "format. Defaults to sdpDir with '_out' appended.");
  basic_options.add_options()(
    "checkpointDir,c", po::value<boost::filesystem::path>(&checkpoint_out),
    "Checkpoints are saved to this directory every checkpointInterval. "
    "Defaults to sdpDir with '.ck' extension.");
  basic_options.add_options()(
    "initialCheckpointDir,i",
    po::value<boost::filesystem::path>(&checkpoint_in),
    "The initial checkpoint directory to load. Defaults to "
    "checkpointDir.");
  basic_options.add_options()(
    "checkpointInterval",
    po::value<int64_t>(&checkpoint_interval)->default_value(3600),
    "Save checkpoints to checkpointDir every checkpointInterval "
    "seconds.");
  basic_options.add_options()(
    "noFinalCheckpoint",
    po::bool_switch(&no_final_checkpoint)->default_value(false),
    "Don't save a final checkpoint after terminating (useful when "
    "debugging).");
  basic_options.add_options()(
    "writeSolution",
    po::value<std::string>(&write_solution_string)->default_value("x,y"s),
    "A comma separated list of vectors and matrices to write into the output "
    "directory.  The default only writes the vectors 'x' and 'y'.  If you add "
    "the 'X' and 'Y' matrices, then the output directory can be used as a "
    "final text "
    "checkpoint.  Runs started from text checkpoints will very close to, but "
    "not bitwise identical to, the original run.\nTo only output the result "
    "(because, for example, you only want to know if SDPB found a primal "
    "feasible point), set this to an empty string.");
  basic_options.add_options()(
    "procGranularity", po::value<size_t>(&proc_granularity)->default_value(1),
    "procGranularity must evenly divide procsPerNode.\n\n"
    "The minimum number of cores in a group, used during load balancing.  "
    "Setting it to anything larger than 1 will make the solution take "
    "longer.  "
    "This option is generally useful only when trying to fit a large problem "
    "in a small machine.");
  basic_options.add_options()("verbosity",
                              po::value<int>(&int_verbosity)->default_value(1),
                              "Verbosity.  0 -> no output, 1 -> regular "
                              "output, 2 -> debug output");

  // We set default parameters using El::BigFloat("1e-10",10)
  // rather than a straight double precision 1e-10 so that results
  // are more reproducible at high precision.  Using double
  // precision defaults results in differences of about 1e-15 in
  // primalObjective after one step.
  po::options_description solver_options("Solver parameters");
  solver_options.add_options()(
    "precision", po::value<size_t>(&precision)->default_value(400),
    "The precision, in the number of bits, for numbers in the "
    "computation. "
    " This should be less than or equal to the precision used when "
    "preprocessing the XML input files with 'pvm2sdp'.  GMP will round "
    "this up to a multiple of 32 or 64, depending on the system.");
  solver_options.add_options()(
    "findPrimalFeasible",
    po::bool_switch(&find_primal_feasible)->default_value(false),
    "Terminate once a primal feasible solution is found.");
  solver_options.add_options()(
    "findDualFeasible",
    po::bool_switch(&find_dual_feasible)->default_value(false),
    "Terminate once a dual feasible solution is found.");
  solver_options.add_options()(
    "detectPrimalFeasibleJump",
    po::bool_switch(&detect_primal_feasible_jump)->default_value(false),
    "Terminate if a primal-step of 1 is taken. This often indicates that "
    "a "
    "primal feasible solution would be found if the precision were high "
    "enough. Try increasing either primalErrorThreshold or precision "
    "and run from the latest checkpoint.");
  solver_options.add_options()(
    "detectDualFeasibleJump",
    po::bool_switch(&detect_dual_feasible_jump)->default_value(false),
    "Terminate if a dual-step of 1 is taken. This often indicates that a "
    "dual feasible solution would be found if the precision were high "
    "enough. Try increasing either dualErrorThreshold or precision "
    "and run from the latest checkpoint.");
  solver_options.add_options()(
    "maxIterations", po::value<int64_t>(&max_iterations)->default_value(500),
    "Maximum number of iterations to run the solver.");
  solver_options.add_options()(
    "maxRuntime",
    po::value<int64_t>(&max_runtime)
      ->default_value(std::numeric_limits<int64_t>::max()),
    "Maximum amount of time to run the solver in seconds.");
  solver_options.add_options()(
    "dualityGapThreshold",
    po::value<El::BigFloat>(&duality_gap_threshold)
      ->default_value(El::BigFloat("1e-30", 10)),
    "Threshold for duality gap (roughly the difference in primal and dual "
    "objective) at which the solution is considered "
    "optimal. Corresponds to SDPA's epsilonStar.");
  solver_options.add_options()(
    "primalErrorThreshold",
    po::value<El::BigFloat>(&primal_error_threshold)
      ->default_value(El::BigFloat("1e-30", 10)),
    "Threshold for feasibility of the primal problem. Corresponds to "
    "SDPA's epsilonBar.");
  solver_options.add_options()(
    "dualErrorThreshold",
    po::value<El::BigFloat>(&dual_error_threshold)
      ->default_value(El::BigFloat("1e-30", 10)),
    "Threshold for feasibility of the dual problem. Corresponds to SDPA's "
    "epsilonBar.");
  solver_options.add_options()(
    "initialMatrixScalePrimal",
    po::value<El::BigFloat>(&initial_matrix_scale_primal)
      ->default_value(El::BigFloat("1e20", 10)),
    "The primal matrix X begins at initialMatrixScalePrimal times the "
    "identity matrix. Corresponds to SDPA's lambdaStar.");
  solver_options.add_options()(
    "initialMatrixScaleDual",
    po::value<El::BigFloat>(&initial_matrix_scale_dual)
      ->default_value(El::BigFloat("1e20", 10)),
    "The dual matrix Y begins at initialMatrixScaleDual times the "
    "identity matrix. Corresponds to SDPA's lambdaStar.");
  solver_options.add_options()(
    "feasibleCenteringParameter",
    po::value<El::BigFloat>(&feasible_centering_parameter)
      ->default_value(El::BigFloat("0.1", 10)),
    "Shrink the complementarity X Y by this factor when the primal and "
    "dual "
    "problems are feasible. Corresponds to SDPA's betaStar.");
  solver_options.add_options()(
    "infeasibleCenteringParameter",
    po::value<El::BigFloat>(&infeasible_centering_parameter)
      ->default_value(El::BigFloat("0.3", 10)),
    "Shrink the complementarity X Y by this factor when either the primal "
    "or dual problems are infeasible. Corresponds to SDPA's betaBar.");
  solver_options.add_options()(
    "stepLengthReduction",
    po::value<El::BigFloat>(&step_length_reduction)
      ->default_value(El::BigFloat("0.7", 10)),
    "Shrink each newton step by this factor (smaller means slower, more "
    "stable convergence). Corresponds to SDPA's gammaStar.");
  solver_options.add_options()(
    "maxComplementarity",
    po::value<El::BigFloat>(&max_complementarity)
      ->default_value(El::BigFloat("1e100", 10)),
    "Terminate if the complementarity mu = Tr(X Y)/dim(X) "
    "exceeds this value.");

  po::options_description cmd_line_options;
  cmd_line_options.add(required_options).add(basic_options).add(solver_options);

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
        }
      else if(variables_map.count("version") != 0)
        {
          if(El::mpi::Rank() == 0)
            {
              std::cout << "SDPB " << SDPB_VERSION_STRING << "\n";
            }
          El::PrintVersion();
          El::PrintConfig();
          El::PrintCCompilerInfo();
          El::PrintCxxCompilerInfo();
        }
      else
        {
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

              po::store(po::parse_config_file(ifs, cmd_line_options),
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

          if(variables_map.count("outDir") == 0)
            {
              out_directory = sdp_directory;
              if(out_directory.filename() == ".")
                {
                  out_directory = out_directory.parent_path();
                }
              out_directory += "_out";
            }

          if(variables_map.count("checkpointDir") == 0)
            {
              checkpoint_out = sdp_directory;
              if(checkpoint_out.filename() == ".")
                {
                  checkpoint_out = checkpoint_out.parent_path();
                }
              checkpoint_out += ".ck";
            }

          if(variables_map.count("initialCheckpointDir") == 0)
            {
              checkpoint_in = checkpoint_out;
            }
          else
            {
              require_initial_checkpoint = true;
            }

          write_solution = Write_Solution(write_solution_string);

          if(El::mpi::Rank() == 0)
            {
              boost::filesystem::create_directories(out_directory);
              boost::filesystem::ofstream ofs(out_directory / "out.txt");
              if(!ofs.good())
                {
                  throw std::runtime_error("Cannot write to outDir: "
                                           + out_directory.string());
                }
            }

          if(int_verbosity != 0 && int_verbosity != 1 && int_verbosity != 2)
            {
              throw std::runtime_error(
                "Invalid number for Verbosity.  Only 0, 1 or 2 are allowed\n");
            }
          else
            {
              verbosity = static_cast<Verbosity>(int_verbosity);
            }
        }
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
