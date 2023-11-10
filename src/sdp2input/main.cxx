#include "../sdp_read.hxx"
#include "../Timers.hxx"

#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

void write_output(const boost::filesystem::path &output_path,
                  Block_File_Format output_format,
                  const std::vector<std::string> &command_arguments,
                  const std::vector<El::BigFloat> &objectives,
                  const std::vector<El::BigFloat> &normalization,
                  const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                  Timers &timers, bool debug);

std::istream &operator>>(std::istream &in, Block_File_Format &format)
{
  std::string token;
  in >> token;
  if(token == "json")
    format = json;
  else if(token == "bin")
    format = bin;
  else
    in.setstate(std::ios_base::failbit);
  return in;
}

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      int precision;
      boost::filesystem::path input_file, output_path;
      Block_File_Format output_format;
      bool debug(false);

      po::options_description options("Basic options");
      options.add_options()("help,h", "Show this helpful message.");
      options.add_options()(
        "input,i", po::value<boost::filesystem::path>(&input_file)->required(),
        "Mathematica, JSON, or NSV file with SDP definition");
      options.add_options()(
        "output,o",
        po::value<boost::filesystem::path>(&output_path)->required(),
        "Directory to place output");
      options.add_options()(
        "precision", po::value<int>(&precision)->required(),
        "The precision, in the number of bits, for numbers in the "
        "computation. ");
      options.add_options()(
        "outputFormat,f",
        po::value<Block_File_Format>(&output_format)
          ->default_value(Block_File_Format::json),
        "Output format for SDP blocks. Could be either 'bin' or 'json'");
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
      Timers timers(debug);
      auto &read_input_timer(timers.add_and_start("read_input"));
      read_input(input_file, objectives, normalization, matrices);
      read_input_timer.stop();
      auto &write_output_timer(timers.add_and_start("write_output"));
      std::vector<std::string> command_arguments;
      for(int arg(0); arg != argc; ++arg)
        {
          command_arguments.emplace_back(argv[arg]);
        }
      write_output(output_path, output_format, command_arguments, objectives,
                   normalization, matrices, timers, debug);
      write_output_timer.stop();
      if(debug)
        {
          timers.write_profile(output_path.string() + ".profiling."
                               + std::to_string(El::mpi::Rank()));
        }
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
