#include "../Approx_Parameters.hxx"

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>

namespace po = boost::program_options;

Approx_Parameters::Approx_Parameters(int argc, char *argv[])
{
  std::string write_solution_string;
  using namespace std::string_literals;

  po::options_description required_options("Required options");
  required_options.add_options()(
    "sdp", po::value<boost::filesystem::path>(&sdp_path)->required(),
    "Directory containing preprocessed SDP data files corresponding to the solution.");
  required_options.add_options()(
    "newSdp", po::value<boost::filesystem::path>(&new_sdp_path)->required(),
    "Directory containing preprocessed data files for the SDP you wish to approximate.");
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
  required_options.add_options()(
    "precision",
    boost::program_options::value<size_t>(&precision)->required(),
    "The precision, in the number of bits, for numbers in the "
    "computation. "
    " This should be less than or equal to the precision used when "
    "preprocessing the XML input files with 'pvm2sdp'.  GMP will round "
    "this up to a multiple of 32 or 64, depending on the system.");

  po::options_description cmd_line_options;
  cmd_line_options.add(required_options);

  po::options_description basic_options("Basic options");
  basic_options.add_options()("help,h", "Show this helpful message.");
  basic_options.add_options()(
    "paramFile,p", po::value<boost::filesystem::path>(&param_path),
    "Any parameter can optionally be set via this file in key=value "
    "format. Command line arguments override values in the parameter "
    "file.");
  basic_options.add_options()(
    "outDir,o", po::value<boost::filesystem::path>(&out_directory),
    "The objective is written this directory in Mathematica "
    "format. Defaults to newSdp with '_out' appended.");
  basic_options.add_options()(
    "procGranularity", po::value<size_t>(&proc_granularity)->default_value(1),
    "procGranularity must evenly divide procsPerNode.\n\n"
    "The minimum number of cores in a group, used during load balancing.  "
    "Setting it to anything larger than 1 will make the solution take "
    "longer.  "
    "This option is generally useful only when trying to fit a large problem "
    "in a small machine.");
  basic_options.add_options()(
                              "solutionDir",
    boost::program_options::value<boost::filesystem::path>(&solution_dir),
    "The directory with the solution of x and y for the primary sdp. Defaults to "
    "primary with '_out' appended.");

  cmd_line_options.add(basic_options);

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
              // TODO: The next line is redundant.  param_file has
              // already been set.  Also, I can use
              // boost::filesystem::ifstream and avoid the
              // .string().c_str() nonsense.
              param_path
                = variables_map["paramFile"].as<boost::filesystem::path>();
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

          if(!boost::filesystem::exists(sdp_path))
            {
              throw std::runtime_error("SDP path '"
                                       + sdp_path.string()
                                       + "' does not exist");
            }

          if(!boost::filesystem::exists(new_sdp_path))
            {
              throw std::runtime_error("New SDP path '"
                                       + new_sdp_path.string()
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

          if(variables_map.count("outDir") == 0)
            {
              out_directory = new_sdp_path;
              if(out_directory.filename() == ".")
                {
                  out_directory = out_directory.parent_path();
                }
              out_directory += "_out";
            }

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
        }
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
