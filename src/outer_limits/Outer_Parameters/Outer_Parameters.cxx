#include "../Outer_Parameters.hxx"

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>

namespace po = boost::program_options;

Outer_Parameters::Outer_Parameters(int argc, char *argv[])
{
  int int_verbosity;
  std::string write_solution_string;
  using namespace std::string_literals;

  po::options_description required_options("Required options");
  required_options.add_options()(
    "functions",
    po::value<boost::filesystem::path>(&functions_path)->required(),
    "Mathematica, JSON, or NSV file with SDP functions evaluated at chebyshev "
    "zeros.");
  required_options.add_options()(
    "points", po::value<boost::filesystem::path>(&points_path)->required(),
    "JSON, or NSV file with initial points.");

  po::options_description cmd_line_options;
  cmd_line_options.add(required_options);

  po::options_description basic_options("Basic options");
  basic_options.add_options()("help,h", "Show this helpful message.");
  basic_options.add_options()("version",
                              "Show version and configuration info.");
  basic_options.add_options()(
    "paramFile,p", po::value<boost::filesystem::path>(&param_path),
    "Any parameter can optionally be set via this file in key=value "
    "format. Command line arguments override values in the parameter "
    "file.");
  basic_options.add_options()(
    "outDir,o", po::value<boost::filesystem::path>(&out_directory),
    "The optimal solution is saved to this directory in Mathematica "
    "format. Defaults to sdp with '_out' appended.");
  basic_options.add_options()("verbosity",
                              po::value<int>(&int_verbosity)->default_value(1),
                              "Verbosity.  0 -> no output, 1 -> regular "
                              "output, 2 -> debug output");

  cmd_line_options.add(basic_options);
  po::options_description solver_options(solver.options());
  solver_options.add_options()(
    "dualityGapReduction",
    po::value<El::BigFloat>(&duality_gap_reduction)
      ->default_value(El::BigFloat("1024", 10)),
    "Shrink the duality gap threshold by this factor during each outer "
    "iteration.  Smaller means slower convergence.");
  cmd_line_options.add(solver_options);

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
              boost::filesystem::ifstream param_file(param_path);
              if(!param_file.good())
                {
                  throw std::runtime_error("Could not open '"
                                           + param_path.string() + "'");
                }

              po::store(po::parse_config_file(param_file, cmd_line_options),
                        variables_map);
            }

          po::notify(variables_map);

          if(!boost::filesystem::exists(functions_path))
            {
              throw std::runtime_error("sdp path '" + functions_path.string()
                                       + "' does not exist");
            }
          if(boost::filesystem::is_directory(functions_path))
            {
              throw std::runtime_error("sdp path '" + functions_path.string()
                                       + "' is a directory, not a file.");
            }

          if(variables_map.count("outDir") == 0)
            {
              out_directory = functions_path.string() + "_out";
            }

          if(variables_map.count("checkpointDir") == 0)
            {
              solver.checkpoint_out = functions_path;
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
