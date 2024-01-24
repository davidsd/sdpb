#include "../Outer_Parameters.hxx"

#include <boost/program_options.hpp>

namespace fs = std::filesystem;
namespace po = boost::program_options;

Outer_Parameters::Outer_Parameters(int argc, char *argv[])
{
  std::string write_solution_string;
  using namespace std::string_literals;

  po::options_description required_options("Required options");
  required_options.add_options()(
    "functions", po::value<fs::path>(&functions_path)->required(),
    "Mathematica, JSON, or NSV file with SDP functions evaluated at Chebyshev "
    "zeros.");
  required_options.add_options()("points",
                                 po::value<fs::path>(&points_path)->required(),
                                 "JSON or NSV file with initial points.");

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
    "out,o", po::value<fs::path>(&output_path),
    "The optimal solution is saved to this file in json "
    "format. Defaults to 'functions' with the ending '_out.json'.");
  basic_options.add_options()("verbosity",
    po::value<Verbosity>(&verbosity)->default_value(Verbosity::regular),
    "Verbosity.  0 -> no output, 1 -> regular output, 2 -> debug output");

  cmd_line_options.add(basic_options);
  po::options_description solver_options(solver.options());
  solver_options.add_options()(
    "dualityGapReduction",
    po::value<El::BigFloat>(&duality_gap_reduction)
      ->default_value(El::BigFloat("1024", 10)),
    "Shrink the duality gap threshold by this factor during each outer "
    "iteration.  Smaller means slower convergence.");
  solver_options.add_options()(
    "meshThreshold",
    po::value<El::BigFloat>(&mesh_threshold)
      ->default_value(El::BigFloat("0.001", 10)),
    "Relative error threshold for when to refine a mesh when approximating a "
    "functional to look for negative regions.");
  solver_options.add_options()(
    "useSVD", po::value<bool>(&use_svd)->default_value(true),
    "Regularize the system using an SVD.  This reduces the required "
    "precision, but can be slow for large problems.");
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
              std::ifstream param_file(param_path);
              if(!param_file.good())
                {
                  throw std::runtime_error("Could not open '"
                                           + param_path.string() + "'");
                }

              po::store(po::parse_config_file(param_file, cmd_line_options),
                        variables_map);
            }

          po::notify(variables_map);

          if(!fs::exists(functions_path))
            {
              throw std::runtime_error("functions path '"
                                       + functions_path.string()
                                       + "' does not exist");
            }
          if(fs::is_directory(functions_path))
            {
              throw std::runtime_error("functions path '"
                                       + functions_path.string()
                                       + "' is a directory, not a file.");
            }

          if(variables_map.count("out") == 0)
            {
              fs::path filename(functions_path);
              filename.replace_extension();
              output_path = filename.string() + "_out.json";
            }

          if(variables_map.count("checkpointDir") == 0)
            {
              solver.checkpoint_out = functions_path;
              if(solver.checkpoint_out.filename() == ".")
                {
                  solver.checkpoint_out = solver.checkpoint_out.parent_path();
                }
              solver.checkpoint_out.replace_extension("ck");
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
              fs::create_directories(output_path.parent_path());
              std::ofstream ofs(output_path);
              if(!ofs.good())
                {
                  throw std::runtime_error("Cannot write to outDir: "
                                           + output_path.string());
                }
            }
        }
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
