#include "../sdp_convert.hxx"

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<boost::filesystem::path> &input_files,
                        boost::filesystem::path &output_dir);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  size_t &num_processed);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  const int rank(El::mpi::Rank());
  try
    {
      int precision;
      std::vector<boost::filesystem::path> input_files;
      boost::filesystem::path output_path;

      // TODO read the following from command line:
      Block_File_Format output_format = json;
      bool debug = false;

      parse_command_line(argc, argv, precision, input_files, output_path);
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
      write_sdpb_input_files(output_path, output_format, rank, num_blocks,
                             command_arguments, objective_const,
                             dual_objective_b, dual_constraint_groups, debug);
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
