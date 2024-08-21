#include "pmp_read/pmp_read.hxx"
#include "pmp2sdp/Block_File_Format.hxx"
#include "pmp2sdp/Dual_Constraint_Group.hxx"
#include "pmp2sdp/write_sdp.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <boost/program_options.hpp>

#include <filesystem>
#include <string>

namespace fs = std::filesystem;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      if(El::mpi::Rank() == 0)
        El::Output("sdp2input is DEPRECATED, please use pmp2sdp instead.");
      // TODO remove sdp2input in 2.8.0 release

      int precision;
      fs::path input_file, output_path;
      Block_File_Format output_format;
      bool debug;

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
      options.add_options()(
        "outputFormat,f",
        po::value<Block_File_Format>(&output_format)
          ->default_value(Block_File_Format::bin),
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
      std::vector<std::string> command_arguments;
      for(int arg(0); arg != argc; ++arg)
        {
          command_arguments.emplace_back(argv[arg]);
        }

      if(variables_map.count("help") != 0)
        {
          std::cout << options << '\n';
          return 0;
        }

      po::notify(variables_map);
      const auto verbosity = debug ? Verbosity::debug : Verbosity::regular;

      // Print command line
      if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
        {
          std::vector<std::string> arg_list(argv, argv + argc);
          for(const auto &arg : arg_list)
            std::cout << arg << " ";
          std::cout << std::endl;
        }

      ASSERT(fs::exists(input_file),
             "Input file does not exist: ", input_file);
      ASSERT(!fs::is_directory(input_file) && input_file != ".",
             "Input file is a directory, not a file:", input_file);

      Environment::set_precision(precision);

      Timers timers(env, verbosity);
      Scoped_Timer timer(timers, "sdp2input");

      auto pmp
        = read_polynomial_matrix_program(env, input_file, verbosity, timers);

      Output_SDP sdp(pmp, command_arguments, timers);
      bool zip = false;
      write_sdp(output_path, sdp, pmp, output_format, zip, timers, verbosity);
      if(El::mpi::Rank() == 0)
        {
          El::Output("Processed ", sdp.num_blocks, " SDP blocks in ",
                     (double)timer.timer().elapsed_milliseconds() / 1000,
                     " seconds, output: ", output_path.string());
        }
      if(verbosity >= Verbosity::debug)
        {
          timers.write_profile(output_path.string() + ".profiling/profiling."
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
