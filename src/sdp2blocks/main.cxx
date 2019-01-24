#include "Positive_Matrix_With_Prefactor.hxx"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

void read_input(const boost::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices);

void write_output_files(
  const boost::filesystem::path &outdir,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  if(El::mpi::Size(El::mpi::COMM_WORLD) != 1)
    {
      if(El::mpi::Rank() == 0)
        {
          std::cerr
            << "sdp2blocks can only be run with a single MPI task, but "
               "was invoked with "
            << El::mpi::Size(El::mpi::COMM_WORLD) << " tasks.\n"
            << std::flush;
        }
      El::Finalize();
      exit(-1);
    }

  try
    {
      int precision;
      boost::filesystem::path input_file, output_dir;

      po::options_description options("Basic options");
      options.add_options()("help,h", "Show this helpful message.");
      options.add_options()(
        "input,i", po::value<boost::filesystem::path>(&input_file)->required(),
        "XML file with SDP definition");
      options.add_options()(
        "output,o",
        po::value<boost::filesystem::path>(&output_dir)->required(),
        "Directory to place output");
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

      if(boost::filesystem::exists(output_dir)
         && !boost::filesystem::is_directory(output_dir))
        {
          throw std::runtime_error("Output directory '" + output_dir.string()
                                   + "' exists and is not a directory");
        }

      mpf_set_default_prec(precision);
      El::gmp::SetPrecision(precision);
      El::mpfr::SetPrecision(precision);
      using Boost_Float = boost::multiprecision::mpfr_float;
      Boost_Float::default_precision(El::gmp::Precision());

      std::vector<El::BigFloat> objectives, normalization;
      std::vector<Positive_Matrix_With_Prefactor> matrices;
      read_input(input_file, objectives, normalization, matrices);
      write_output_files(output_dir, objectives, normalization, matrices);
    }
  catch(std::runtime_error &e)
    {
      std::cerr << "Error: " << e.what() << "\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
  catch(po::error &e)
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
