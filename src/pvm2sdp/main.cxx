#include "sdp_convert/sdp_convert.hxx"

namespace fs = std::filesystem;

// TODO use boost::program_options
void parse_command_line(int argc, char **argv,
                        Block_File_Format &output_format, int &precision,
                        std::vector<fs::path> &input_files,
                        fs::path &output_dir);

void read_input_files(
  const std::vector<fs::path> &input_files, El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  size_t &num_processed);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  const int rank(El::mpi::Rank());
  try
    {
      Block_File_Format output_format = bin;
      int precision;
      std::vector<fs::path> input_files;
      fs::path output_path;

      bool debug = false;

      parse_command_line(argc, argv, output_format, precision, input_files,
                         output_path);
      El::gmp::SetPrecision(precision);

      El::BigFloat objective_const;
      std::vector<El::BigFloat> dual_objective_b;
      std::vector<Dual_Constraint_Group> dual_constraint_groups;
      size_t num_blocks(0);
      read_input_files(input_files, objective_const, dual_objective_b,
                       dual_constraint_groups, num_blocks);

      std::vector<std::string> command_arguments;
      for(int arg(0); arg != argc; ++arg)
        {
          command_arguments.emplace_back(argv[arg]);
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
      Timers timers(env, debug);
      write_sdpb_input_files(output_path, output_format, rank, num_blocks,
                             command_arguments, objective_const,
                             dual_objective_b, dual_constraint_groups, timers,
                             debug);
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
