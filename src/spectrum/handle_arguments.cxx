#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/Environment.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <boost/program_options.hpp>
#include <filesystem>
#include <optional>

namespace fs = std::filesystem;

void handle_arguments(const int &argc, char **argv,
                      std::optional<Boost_Float> &threshold,
                      El::BigFloat &max_zero, El::BigFloat &min_zero_distance,
                      fs::path &pmp_info_path, fs::path &solution_dir,
                      fs::path &c_minus_By_path, fs::path &output_path,
                      bool &need_lambda,
                      std::optional<El::BigFloat> &min_eigenvalue_ratio,
                      Verbosity &verbosity)
{
  int precision;
  std::string threshold_string, max_zero_string, min_zero_distance_string,
    mesh_threshold_string, format_string, min_eigenvalue_ratio_string;

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
    "threshold", po::value<std::string>(&threshold_string),
    "Threshold for when a functional is considered to be zero.\n"
    "By default, --threshold=sqrt(dualityGap).");
  options.add_options()(
    "output,o", po::value<fs::path>(&output_path)->required(), "Output file");
  options.add_options()(
    "precision", po::value<int>(&precision)->required(),
    "The precision, in the number of bits, for numbers in the "
    "computation. ");
  options.add_options()(
    "maxZero,m", po::value<std::string>(&max_zero_string)->default_value("0"),
    "Spectrum will ignore all zeros larger than --maxZero. "
    "--maxZero=0 means no limit.");
  options.add_options()(
    // TODO: set default distance finite value, e.g. 2^{-precision/2}?
    "minZeroDistance",
    po::value<std::string>(&min_zero_distance_string)->default_value("0"),
    "Several zeros within --minZeroDistance from each other "
    "are replaced with their arithmetic mean. ");
  options.add_options()("lambda",
                        po::value<bool>(&need_lambda)->default_value(true),
                        "If true, compute Λ and its associated error.");
  options.add_options()(
    "minEigenvalueRatio", po::value<std::string>(&min_eigenvalue_ratio_string),
    "When computing Λ, keep only eigenvalues larger than "
    "minEigenvalueRatio * max(eigenvalues).\n"
    "To filter out numerical noise, set this value "
    "somewhat higher than dualityGap of your SDPB solution.\n"
    "By default, --minEigenvalueRatio=sqrt(dualityGap).");
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
  positional.add("maxZero", 1);

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

  // Set parameters
  {
    Environment::set_precision(precision);
    if(variables_map.count("threshold") != 0)
      threshold.emplace(Boost_Float(threshold_string));
    max_zero = El::BigFloat(max_zero_string);
    min_zero_distance = El::BigFloat(min_zero_distance_string);
    if(c_minus_By_path.empty())
      c_minus_By_path = solution_dir / "c_minus_By" / "c_minus_By.json";
    if(variables_map.count("minEigenvalueRatio") != 0)
      min_eigenvalue_ratio.emplace(min_eigenvalue_ratio_string);
  }

  // Asserts and warnings
  if(El::mpi::Rank() == 0)
    {
      ASSERT(max_zero >= 0, "--maxZero must be non-negative");
      ASSERT(min_zero_distance >= 0, "--minZeroDistance must be non-negative");

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

      if(variables_map.count("solution") == 0)
        {
          ASSERT(need_lambda == false,
                 "--solution must be specified unless --lambda=false");
          ASSERT(variables_map.count("cMinusBy") != 0,
                 "Please specify either --solution or --cMinusBy");
          ASSERT(threshold.has_value(),
                 "To set --threshold automatically, "
                 "you need to specify --solution directory");
        }

      if(variables_map.count("cMinusBy") == 0 || need_lambda)
        {
          ASSERT(fs::exists(solution_dir),
                 "--solution directory does not exist: ", solution_dir);
          ASSERT(fs::is_directory(solution_dir),
                 "--solution is not a directory: ", solution_dir);
        }

      if(min_eigenvalue_ratio.has_value())
        {
          if(!need_lambda)
            {
              PRINT_WARNING(
                "--minEigenvalueRatio will be ignored since --lambda=false");
            }

          ASSERT(min_eigenvalue_ratio.value() >= 0
                   && min_eigenvalue_ratio.value() <= 1,
                 "--minEigenvalueRatio=", min_eigenvalue_ratio_string,
                 " should be in range [0,1].");
        }

      ASSERT(fs::exists(c_minus_By_path), DEBUG_STRING(c_minus_By_path));
      ASSERT(fs::is_regular_file(c_minus_By_path),
             DEBUG_STRING(c_minus_By_path));

      ASSERT(
        // pmp_info.json is a regular file in sdp directory:
        fs::exists(pmp_info_path)
          // or pmp_info.json is inside sdp.zip archive:
          || (fs::exists(pmp_info_path.parent_path())
              && !fs::is_directory(pmp_info_path.parent_path())),
        "--pmpInfo file does not exist: ", pmp_info_path);
      ASSERT(!fs::is_directory(pmp_info_path),
             "--pmpInfo path is a directory, not a file: ", pmp_info_path);

      ASSERT(output_path != ".", "Output file is a directory: ", output_path);
      ASSERT(!(fs::exists(output_path) && fs::is_directory(output_path)),
             "Output file exists and is a directory: ", output_path);
    }

  // Wait for checks performed at rank=0
  El::mpi::Barrier();
}
