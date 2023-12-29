#include "../Approx_Parameters.hxx"

#include <boost/program_options.hpp>

namespace fs = std::filesystem;
namespace po = boost::program_options;

Approx_Parameters::Approx_Parameters(int argc, char *argv[])
{
  std::string write_solution_string;
  using namespace std::string_literals;

  po::options_description required_options("Required options");
  required_options.add_options()(
    "sdp", po::value<fs::path>(&sdp_path)->required(),
    "File or directory containing preprocessed SDP input corresponding to the "
    "solution.");
  required_options.add_options()(
    "precision", boost::program_options::value<size_t>(&precision)->required(),
    "The precision, in the number of bits, for numbers in the "
    "computation. "
    " This should be less than or equal to the precision used when "
    "generating the solution with 'sdpb'.  GMP will round "
    "this up to a multiple of 32 or 64, depending on the system.");

  po::options_description cmd_line_options;
  cmd_line_options.add(required_options);

  po::options_description basic_options("Basic options");
  basic_options.add_options()("help,h", "Show this helpful message.");
  basic_options.add_options()(
    "paramFile,p", po::value<fs::path>(&param_path),
    "Any parameter can optionally be set via this file in key=value "
    "format. Command line arguments override values in the parameter "
    "file.");
  basic_options.add_options()(
    "newSdp", po::value<fs::path>(&new_sdp_path),
    "A file containing either preprocessed input for the SDP you wish to "
    "approximate, or a null separated list of files with preprocessed input.");
  basic_options.add_options()(
    "procGranularity", po::value<size_t>(&proc_granularity)->default_value(1),
    "procGranularity must evenly divide procsPerNode.\n\n"
    "The minimum number of cores in a group, used during load balancing.  "
    "Setting it to anything larger than 1 will make the solution take "
    "longer.  "
    "This option is generally useful only when trying to fit a large problem "
    "in a small machine.");
  basic_options.add_options()(
    "solutionDir", boost::program_options::value<fs::path>(&solution_dir),
    "The directory with the text format solutions of x and y for the primary "
    "sdp. It must also contain either X and Y or the solver state. "
    "Defaults to sdp with '_out' appended.");
  basic_options.add_options()(
    "writeSolverState",
    po::bool_switch(&write_solver_state)->default_value(false),
    "Write the solver state in solutionDir.  This allows later invocations of "
    "approx_objective to skip the time consuming part of setting up the "
    "solver.");
  basic_options.add_options()(
    "linear",
    po::bool_switch(&linear_only)->default_value(false),
    "Only compute the linear correction, not the quadratic correction.  "
    "This avoids having to compute an expensive inverse.");

  cmd_line_options.add(basic_options);

  po::options_description obsolete_options("Obsolete options");
  obsolete_options.add_options()(
    "procsPerNode", po::value<size_t>(),
    "[OBSOLETE] The number of MPI processes running on a node. "
    "Determined automatically from MPI environment.");
  cmd_line_options.add(obsolete_options);

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

          if(variables_map.count("procsPerNode") != 0)
            {
              El::Output("--procsPerNode option is obsolete. The number of "
                         "MPI processes running on a node is determined "
                         "automatically from MPI environment.");
            }

          if(!fs::exists(sdp_path))
            {
              throw std::runtime_error("SDP path '" + sdp_path.string()
                                       + "' does not exist");
            }

          if(!new_sdp_path.empty() && !fs::exists(new_sdp_path))
            {
              throw std::runtime_error("New SDP path '" + new_sdp_path.string()
                                       + "' does not exist");
            }

          if(variables_map.count("solutionDir") == 0)
            {
              solution_dir = sdp_path;
              if(solution_dir.filename() == ".")
                {
                  solution_dir = solution_dir.parent_path();
                }
              solution_dir += "_out";
            }
        }
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
