#include "../Dynamical_Parameters.hxx"

#include <El.hpp>

#include <boost/program_options.hpp>
#include <fstream>

namespace po = boost::program_options;

Dynamical_Parameters::Dynamical_Parameters(int argc, char *argv[])
{
  int int_verbosity;
  std::string write_solution_string;
  using namespace std::string_literals;

  po::options_description required_options("Required options");
  required_options.add_options()(
    "sdpDir,s", po::value<std::filesystem::path>(&sdp_path)->required(),
    "Directory containing the preprocessed centering SDP data files.");
  //required_options.add_options()(
  //  "newSdpDirs", po::value<std::filesystem::path>(&new_sdp_path)->required(),
  //  "Directory containing the preprocessed SDP data files around the center SDP in external parameter space.");

  po::options_description cmd_line_options;
  cmd_line_options.add(required_options);

  po::options_description basic_options("Basic options");
  basic_options.add_options()("help,h", "Show this helpful message.");
  basic_options.add_options()("version",
                              "Show version and configuration info.");
  basic_options.add_options()(
    "paramFile,p", po::value<std::filesystem::path>(&param_path),
    "Any parameter can optionally be set via this file in key=value "
    "format. Command line arguments override values in the parameter "
    "file.");
  basic_options.add_options()(
    "outDir,o", po::value<std::filesystem::path>(&out_directory),
    "The optimal solution is saved to this directory in Mathematica "
    "format. Defaults to sdpDir with '_out' appended.");
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
    "procGranularity must evenly divide number of processes per node.\n\n"
    "The minimum number of cores in a group, used during load balancing.  "
    "Setting it to anything larger than 1 will make the solution take "
    "longer.  "
    "This option is generally useful only when trying to fit a large problem "
    "in a small machine.");
  basic_options.add_options()("verbosity",
                              po::value<int>(&int_verbosity)->default_value(1),
                              "Verbosity.  0 -> no output, 1 -> regular "
                              "output, 2 -> debug output");

  po::options_description obsolete_options("Obsolete options");
  obsolete_options.add_options()(
    "procsPerNode", po::value<size_t>(),
    "[OBSOLETE] The number of MPI processes running on a node. "
    "Determined automatically from MPI environment.");

  cmd_line_options.add(basic_options);
  cmd_line_options.add(obsolete_options);
  cmd_line_options.add(solver.options());

  po::variables_map variables_map;
  try
    {
      po::store(
        po::parse_command_line(argc, argv, cmd_line_options,
                               po::command_line_style::unix_style
                                 ^ po::command_line_style::allow_short),
        variables_map);
      //po::store(po::parse_command_line(argc, argv, cmd_line_options),
      //  variables_map);

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
              // Where is SDPB_VERSION_STRING initialized?
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
              // TODO: The next line is redundant.  param_file has
              // already been set.  Also, I can use
              // std::ifstream and avoid the
              // .string().c_str() nonsense.
              param_path
                = variables_map["paramFile"].as<std::filesystem::path>();
              std::ifstream ifs(param_path.string().c_str());
              if(!ifs.good())
                {
                  throw std::runtime_error("Could not open '"
                                           + param_path.string() + "'");
                }

              po::store(po::parse_config_file(ifs, cmd_line_options),
                        variables_map);
            }

          po::notify(variables_map);

          if(!std::filesystem::exists(sdp_path))
            {
              throw std::runtime_error("sdp directory '" + sdp_path.string()
                                       + "' does not exist");
            }
          if(!std::filesystem::exists(solver.new_sdp_path))
            {
              throw std::runtime_error("new sdp directory '"
                                       + solver.new_sdp_path.string()
                                       + "' does not exist");
            }
          if(variables_map.count("outDir") == 0)
            {
              out_directory = sdp_path;
              if(out_directory.filename() == ".")
                {
                  out_directory = out_directory.parent_path();
                }
              out_directory += "_out";
            }

          if(variables_map.count("checkpointDir") == 0)
            {
              solver.solver_parameters.checkpoint_out = sdp_path;
              if(solver.solver_parameters.checkpoint_out.filename() == ".")
                {
                  solver.solver_parameters.checkpoint_out
                    = solver.solver_parameters.checkpoint_out.parent_path();
                }
              solver.solver_parameters.checkpoint_out += ".ck";
            }

          if(variables_map.count("initialCheckpointDir") == 0)
            {
              solver.solver_parameters.checkpoint_in
                = solver.solver_parameters.checkpoint_out;
            }
          else
            {
              require_initial_checkpoint = true;
            }

          write_solution = Write_Solution(write_solution_string);

          if(El::mpi::Rank() == 0)
            {
              std::filesystem::create_directories(out_directory);
              std::ofstream ofs(out_directory / "out.txt");
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

          if(solver.n_external_parameters == 0)
            {
              throw std::runtime_error(
                "The number of external parameters to be varied is zero.");
            }

          if(variables_map.count("oldSchurDir") == 0)
            {
              solver.external_corrector_Q = false;
            }

          if(variables_map.count("prevGradientBFGS") == 0
             || variables_map.count("prevExternalStep") == 0)
            solver.prev_grad_step_validQ = false;
          else
            solver.prev_grad_step_validQ = true;

          if(variables_map.count("prevHessianBFGS") == 0)
            solver.use_Hmixed_for_BFGS = true;
          else
            solver.use_Hmixed_for_BFGS = false;

          if(variables_map.count("findBoundaryDirection") == 0)
            solver.find_boundary = false;
          else
            solver.find_boundary = true;

          if(solver.navigatorWithLogDetX)
            solver.beta_for_mu_logdetX = El::BigFloat("-1");
          else
            solver.beta_for_mu_logdetX = El::BigFloat("0");
        }
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
