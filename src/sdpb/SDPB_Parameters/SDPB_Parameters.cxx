#include "../SDPB_Parameters.hxx"

#include <boost/program_options.hpp>

namespace fs = std::filesystem;

namespace po = boost::program_options;

SDPB_Parameters::SDPB_Parameters(int argc, char *argv[])
{
  std::string write_solution_string;
  using namespace std::string_literals;

  po::options_description required_options("Required options");
  required_options.add_options()(
    "sdpDir,s", po::value<fs::path>(&sdp_path)->required(),
    "Directory containing preprocessed SDP data files.");

  po::options_description cmd_line_options;
  cmd_line_options.add(required_options);

  po::options_description basic_options("Basic options");
  basic_options.add_options()("help,h", "Show this helpful message.");
  basic_options.add_options()("version",
                              "Show version and configuration info.");
  basic_options.add_options()(
    "paramFile,p", po::value<fs::path>(&param_path),
    "Any parameter can optionally be set via this file in key=value "
    "format. Command line arguments override values in the parameter "
    "file.");
  basic_options.add_options()(
    "outDir,o", po::value<fs::path>(&out_directory),
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
    po::value<Verbosity>(&verbosity)->default_value(Verbosity::regular),
    "Verbosity.  0 -> no output, 1 -> regular output, 2 -> debug output");

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
      po::store(po::parse_command_line(argc, argv, cmd_line_options),
                variables_map);

      if(variables_map.count("help") != 0)
        {
          if(El::mpi::Rank() == 0)
            {
		    std::cout << "SDPB v" << SDPB_VERSION_STRING << "\n";
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
              param_path = variables_map["paramFile"].as<fs::path>();
              std::ifstream ifs(param_path);
              if(!ifs.good())
                {
                  throw std::runtime_error("Could not open '"
                                           + param_path.string() + "'");
                }

              po::store(po::parse_config_file(ifs, cmd_line_options),
                        variables_map);
            }

          po::notify(variables_map);

          if(!fs::exists(sdp_path))
            {
              throw std::runtime_error("sdp directory '"
                                       + sdp_path.string()
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
              solver.checkpoint_out = sdp_path;
              if(solver.checkpoint_out.filename() == ".")
                {
                  solver.checkpoint_out = solver.checkpoint_out.parent_path();
                }
              solver.checkpoint_out += ".ck";
            }

          if(variables_map.count("initialCheckpointDir") == 0)
            {
              solver.checkpoint_in = solver.checkpoint_out;
            }
          else
            {
              require_initial_checkpoint = true;
            }

          write_solution = Write_Solution(write_solution_string);

          if(El::mpi::Rank() == 0)
            {
              fs::create_directories(out_directory);
              std::ofstream ofs(out_directory / "out.txt");
              if(!ofs.good())
                {
                  throw std::runtime_error("Cannot write to outDir: "
                                           + out_directory.string());
                }
            }

          if(El::mpi::Rank() == 0 && verbosity >= Verbosity::regular
             && variables_map.count("procsPerNode") != 0)
            {
              El::Output("--procsPerNode option is obsolete. The number of "
                         "MPI processes running on a node is determined "
                         "automatically from MPI environment.");
            }
        }
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
