#include "pmp_read/pmp_read.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/program_options.hpp>
#include <filesystem>

namespace fs = std::filesystem;
namespace po = boost::program_options;

void write_functions(
  const fs::path &output_path, const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<Polynomial_Vector_Matrix> &matrices);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      // TODO fix parallel
      if(El::mpi::Size() > 1)
        RUNTIME_ERROR("sdp2functions cannot work in parallel!");

      int precision;
      fs::path input_file, output_path;
      bool debug(false);

      po::options_description options("Basic options");
      options.add_options()("help,h", "Show this helpful message.");
      options.add_options()(
        "input,i", po::value<fs::path>(&input_file)->required(),
        "Mathematica, JSON, or NSV file with SDP definition");
      options.add_options()("output,o",
                            po::value<fs::path>(&output_path)->required(),
                            "Directory to place output");
      options.add_options()(
        "precision", po::value<int>(&precision)->required(),
        "The precision, in the number of bits, for numbers in the "
        "computation. ");
      options.add_options()("debug",
                            po::value<bool>(&debug)->default_value(false),
                            "Write out debugging output.");

      po::positional_options_description positional;
      positional.add("precision", 1);
      positional.add("input", 1);
      positional.add("output", 1);

      po::variables_map variables_map;
      po::store(po::parse_command_line(argc, argv, options), variables_map);

      if(variables_map.count("help") != 0)
        {
          std::cout << options << '\n';
          return 0;
        }

      po::notify(variables_map);

      ASSERT(fs::exists(input_file),
             "Input file does not exist: ", input_file);
      ASSERT(!fs::is_directory(input_file) && input_file != ".",
             "Input file is a directory, not a file:", input_file);
      ASSERT(output_path != ".", "Output file is a directory: ", output_path);
      ASSERT(!(fs::exists(output_path) && fs::is_directory(output_path)),
             "Output file exists and is a directory: ", output_path);

      El::gmp::SetPrecision(precision);
      // El::gmp wants base-2 bits, but boost::multiprecision wants
      // base-10 digits.
      Boost_Float::default_precision(precision * log(2) / log(10));

      const auto pmp = read_polynomial_matrix_program(input_file);
      write_functions(output_path, pmp.objective, pmp.normalization,
                      pmp.matrices);
    }
  catch(std::exception &e)
    {
      std::cerr << "Error: " << e.what() << "\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
  catch(...)
    {
      std::cerr << "Unknown Error\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
