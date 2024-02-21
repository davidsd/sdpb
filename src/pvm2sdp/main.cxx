#include "pmp_read/pmp_read.hxx"
#include "pmp2sdp/write_sdp.hxx"

namespace fs = std::filesystem;

// TODO use boost::program_options
void parse_command_line(int argc, char **argv,
                        Block_File_Format &output_format, int &precision,
                        std::vector<fs::path> &input_files,
                        fs::path &output_dir,
                        std::vector<std::string> &command_arguments);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      if(El::mpi::Rank() == 0)
        El::Output("pvm2sdp is DEPRECATED, please use pmp2sdp instead.");
      // TODO remove pvm2sdp in 2.8.0 release

      Block_File_Format output_format = bin;
      int precision;
      std::vector<fs::path> input_files;
      fs::path output_path;
      std::vector<std::string> command_arguments;

      bool debug = false;

      parse_command_line(argc, argv, output_format, precision, input_files,
                         output_path, command_arguments);
      El::gmp::SetPrecision(precision);

      Timers timers(env, debug);

      auto pmp = read_polynomial_matrix_program(env, input_files, timers);
      Output_SDP sdp(pmp, command_arguments, timers);
      bool zip = false;
      write_sdp(output_path, sdp, output_format, zip, timers, debug);
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
