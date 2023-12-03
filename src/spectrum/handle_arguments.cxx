#include "Format.hxx"
#include "sdpb_util/Boost_Float.hxx"

#include <El.hpp>

#include <boost/program_options.hpp>
#include <filesystem>

namespace fs = std::filesystem;

void handle_arguments(const int &argc, char **argv, El::BigFloat &threshold,
                      El::BigFloat &mesh_threshold, Format &format,
                      fs::path &input_path, fs::path &solution_dir,
                      fs::path &output_path, bool &need_lambda)
{
  int precision;
  std::string threshold_string, mesh_threshold_string, format_string;

  namespace po = boost::program_options;

  po::options_description options("Basic options");
  options.add_options()("help,h", "Show this helpful message.");
  options.add_options()("input,i",
                        po::value<fs::path>(&input_path)->required(),
                        "Mathematica, JSON, or NSV file with SDP definition");
  options.add_options()(
    "solution", po::value<fs::path>(&solution_dir)->required(),
    "SDPB output directory containing the solutions for y and x (e.g. "
    "'y.txt')");
  options.add_options()(
    "threshold", po::value<std::string>(&threshold_string)->required(),
    "Threshold for when a functional is considered to be zero.");
  options.add_options()(
    "meshThreshold",
    po::value<std::string>(&mesh_threshold_string)->default_value("0.001"),
    "Relative error threshold for when to refine a mesh when approximating a "
    "functional to look for zeros.");
  options.add_options()("format",
                        po::value<std::string>(&format_string)->required(),
                        "Format of input file: Either PVM (Polynomial Vector "
                        "Matrix), or PMP (Positive Matrix with Prefactor).");
  options.add_options()(
    "output,o", po::value<fs::path>(&output_path)->required(), "Output file");
  options.add_options()(
    "precision", po::value<int>(&precision)->required(),
    "The precision, in the number of bits, for numbers in the "
    "computation. ");
  options.add_options()("lambda",
                        po::value<bool>(&need_lambda)->default_value(true),
                        "If true, compute Λ and its associated error.");

  po::positional_options_description positional;
  positional.add("precision", 1);
  positional.add("input", 1);
  positional.add("output", 1);

  po::variables_map variables_map;
  po::store(po::parse_command_line(argc, argv, options), variables_map);

  if(variables_map.count("help") != 0)
    {
      std::cout << options << '\n';
      exit(0);
    }

  po::notify(variables_map);

  if(!fs::exists(input_path))
    {
      throw std::runtime_error("Input file '" + input_path.string()
                               + "' does not exist");
    }
  if(fs::is_directory(input_path))
    {
      throw std::runtime_error("Input file '" + input_path.string()
                               + "' is a directory, not a file");
    }
  if(!fs::exists(solution_dir))
    {
      throw std::runtime_error("Solution file '" + solution_dir.string()
                               + "' does not exist");
    }

  if(output_path == ".")
    {
      throw std::runtime_error("Output file '" + output_path.string()
                               + "' is a directory");
    }
  if(fs::exists(output_path) && fs::is_directory(output_path))
    {
      throw std::runtime_error("Output file '" + output_path.string()
                               + "' exists and is a directory");
    }

  El::gmp::SetPrecision(precision);
  // El::gmp wants base-2 bits, but boost::multiprecision wants
  // base-10 digits.
  Boost_Float::default_precision(precision * log(2) / log(10));

  threshold = El::BigFloat(threshold_string);
  mesh_threshold = El::BigFloat(mesh_threshold_string);

  if(format_string == "PVM")
    {
      format = Format::Polynomial_Vector_Matrix;
    }
  else if(format_string == "PMP")
    {
      format = Format::Positive_Matrix_with_Prefactor;
    }
  else
    {
      throw std::runtime_error(
        "Unknown format.  Expected 'PVM' or 'PMP', but found '" + format_string
        + "'");
    }
}
