#include "../sdp_read.hxx"
#include "../sdp_solve.hxx"
#include "../read_vector.hxx"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

void write_spectrum(const boost::filesystem::path &output_path,
                    const std::vector<std::vector<El::BigFloat>> &zeros);

namespace po = boost::program_options;

std::vector<std::vector<El::BigFloat>>
compute_spectrum(const std::vector<El::BigFloat> &normalization,
                 const El::DistMatrix<El::BigFloat> &y,
                 const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                 const El::BigFloat &threshold);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      int precision;
      std::string threshold_string;
      boost::filesystem::path input_file, solution_file, output_path;

      po::options_description options("Basic options");
      options.add_options()("help,h", "Show this helpful message.");
      options.add_options()(
        "input,i", po::value<boost::filesystem::path>(&input_file)->required(),
        "Mathematica, JSON, or NSV file with SDP definition");
      options.add_options()(
        "solution",
        po::value<boost::filesystem::path>(&solution_file)->required(),
        "SDPB output file containing the solution for y (e.g. 'y.txt')");
      options.add_options()(
        "threshold", po::value<std::string>(&threshold_string)->required(),
        "Threshold for when a functional is considered to be zero.");
      options.add_options()(
        "output,o",
        po::value<boost::filesystem::path>(&output_path)->required(),
        "Output file");
      options.add_options()(
        "precision", po::value<int>(&precision)->required(),
        "The precision, in the number of bits, for numbers in the "
        "computation. ");

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

      if(!boost::filesystem::exists(input_file))
        {
          throw std::runtime_error("Input file '" + input_file.string()
                                   + "' does not exist");
        }
      if(boost::filesystem::is_directory(input_file))
        {
          throw std::runtime_error("Input file '" + input_file.string()
                                   + "' is a directory, not a file");
        }
      if(!boost::filesystem::exists(solution_file))
        {
          throw std::runtime_error("Solution file '" + solution_file.string()
                                   + "' does not exist");
        }
      if(boost::filesystem::is_directory(solution_file))
        {
          throw std::runtime_error("Solution file '" + solution_file.string()
                                   + "' is a directory, not a file");
        }

      if(output_path.filename_is_dot())
        {
          throw std::runtime_error("Output file '" + output_path.string()
                                   + "' is a directory");
        }
      if(boost::filesystem::exists(output_path)
         && boost::filesystem::is_directory(output_path))
        {
          throw std::runtime_error("Output file '" + output_path.string()
                                   + "' exists and is a directory");
        }

      El::gmp::SetPrecision(precision);
      // El::gmp wants base-2 bits, but boost::multiprecision wants
      // base-10 digits.
      Boost_Float::default_precision(precision * log(2) / log(10));

      std::vector<El::BigFloat> objectives, normalization;
      std::vector<Positive_Matrix_With_Prefactor> matrices;
      read_input(input_file, objectives, normalization, matrices);
      El::DistMatrix<El::BigFloat> y(normalization.size() - 1, 1);
      read_text_block(y, solution_file);
      const std::vector<std::vector<El::BigFloat>> zeros(compute_spectrum(
        normalization, y, matrices, El::BigFloat(threshold_string)));
      write_spectrum(output_path, zeros);
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
