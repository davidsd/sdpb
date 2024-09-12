#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/Environment.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <boost/program_options.hpp>
#include <filesystem>

namespace fs = std::filesystem;

void handle_arguments(const int &argc, char **argv, El::BigFloat &threshold,
                      fs::path &pmp_info_path, fs::path &solution_dir,
                      fs::path &c_minus_By_path, fs::path &output_path,
                      bool &need_lambda, Verbosity &verbosity)
{
  int precision;
  std::string threshold_string, mesh_threshold_string, format_string;

  namespace po = boost::program_options;

  po::options_description options("Basic options");
  options.add_options()("help,h", "Show this helpful message.");
  options.add_options()(
    "pmpInfo,i", po::value<fs::path>(&pmp_info_path)->required(),
    "pmp_info.json with relevant information about PMP blocks. "
    "This file is written to SDP directory by pmp2sdp.");
  options.add_options()(
    "solution", po::value<fs::path>(&solution_dir),
    "SDPB output directory containing the vectors c_minus_By "
    "(file 'c_minus_By/c_minus_By.json') and x (files 'x_*.txt').\n"
    "If --lambda=false, you may omit --solution and specify only --cMinusBy.");
  options.add_options()(
    "cMinusBy", po::value<fs::path>(&c_minus_By_path),
    "Path to c_minus_By.json with the block vector (c - B.y). "
    "By default, equals to '${--solution}/c_minus_By/c_minus_By.json'.");
  options.add_options()(
    "threshold", po::value<std::string>(&threshold_string)->required(),
    "Threshold for when a functional is considered to be zero.");
  options.add_options()(
    "output,o", po::value<fs::path>(&output_path)->required(), "Output file");
  options.add_options()(
    "precision", po::value<int>(&precision)->required(),
    "The precision, in the number of bits, for numbers in the "
    "computation. ");
  options.add_options()("lambda",
                        po::value<bool>(&need_lambda)->default_value(true),
                        "If true, compute Λ and its associated error.");
  options.add_options()(
    "verbosity",
    po::value<Verbosity>(&verbosity)->default_value(Verbosity::regular),
    "Verbosity.  0 -> no output, 1 -> regular output, 2 -> debug output, 3 -> "
    "trace output");

  options.add_options()(
    "meshThreshold", po::value<std::string>(&mesh_threshold_string),
    "[OBSOLETE] Relative error threshold for when to refine a mesh when "
    "approximating a functional to look for zeros.");
  options.add_options()(
    "format", po::value<std::string>(&format_string),
    "[OBSOLETE] Format of input file. Determined automatically.");

  po::positional_options_description positional;
  positional.add("precision", 1);
  positional.add("pmpInfo", 1);
  positional.add("solution", 1);
  positional.add("output", 1);
  positional.add("threshold", 1);

  po::variables_map variables_map;
  po::store(po::command_line_parser(argc, argv)
              .options(options)
              .positional(positional)
              .run(),
            variables_map);

  if(variables_map.count("help") != 0)
    {
      std::cout << options << '\n';
      exit(0);
    }
  po::notify(variables_map);

  // Print warnings
  if(El::mpi::Rank() == 0)
    {
      if(variables_map.count("format") != 0)
        {
          PRINT_WARNING("--format option is obsolete. Input file format is "
                        "determined automatically.");
        }
      if(variables_map.count("meshThreshold") != 0)
        {
          PRINT_WARNING(
            "--meshThreshold option is obsolete and will be ignored");
        }
    }
  if(variables_map.count("solution") == 0)
    {
      ASSERT(variables_map.count("cMinusBy") != 0);
    }

  if(c_minus_By_path.empty() || need_lambda)
    {
      ASSERT(fs::exists(solution_dir),
             "--solution directory does not exist: ", solution_dir);
      ASSERT(fs::is_directory(solution_dir),
             "--solution is not a directory: ", solution_dir);
    }
  if(c_minus_By_path.empty())
    {
      c_minus_By_path = solution_dir / "c_minus_By" / "c_minus_By.json";
    }

  ASSERT(fs::exists(c_minus_By_path), DEBUG_STRING(c_minus_By_path));
  ASSERT(fs::is_regular_file(c_minus_By_path), DEBUG_STRING(c_minus_By_path));

  ASSERT(
    // pmp_info.json is a regular file in sdp directory:
    fs::exists(pmp_info_path)
      // or pmp_info.json is inside sdp.zip archive:
      || (fs::exists(pmp_info_path.parent_path())
          && !fs::is_directory(pmp_info_path.parent_path())),
    "--pmpInfo file does not exist: ", pmp_info_path);
  ASSERT(!fs::is_directory(pmp_info_path),
         "--pmpInfo path is a directory, not a file: ", pmp_info_path);

  ASSERT(fs::exists(solution_dir),
         "Solution directory does not exist: ", solution_dir);

  ASSERT(output_path != ".", "Output file is a directory: ", output_path);
  ASSERT(!(fs::exists(output_path) && fs::is_directory(output_path)),
         "Output file exists and is a directory: ", output_path);

  Environment::set_precision(precision);

  threshold = El::BigFloat(threshold_string);
}
