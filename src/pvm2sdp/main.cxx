#include "../sdp_convert.hxx"

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<boost::filesystem::path> &input_files,
                        boost::filesystem::path &output_dir);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  std::vector<size_t> &indices, size_t &num_processed);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  const int rank(El::mpi::Rank());
  try
    {
      int precision;
      std::vector<boost::filesystem::path> input_files;
      boost::filesystem::path output_dir;

      parse_command_line(argc, argv, precision, input_files, output_dir);
      El::gmp::SetPrecision(precision);

      std::vector<size_t> indices;
      El::BigFloat objective_const;
      std::vector<El::BigFloat> dual_objective_b;
      std::vector<Dual_Constraint_Group> dual_constraint_groups;
      size_t num_blocks(0);
      read_input_files(input_files, objective_const, dual_objective_b,
                       dual_constraint_groups, indices, num_blocks);

      std::vector<std::string> command_arguments;
      for(int arg(0); arg != argc; ++arg)
        {
          command_arguments.emplace_back(argv[arg]);
        }
      write_sdpb_input_files(output_dir, rank, num_blocks,
                             command_arguments, indices, objective_const,
                             dual_objective_b, dual_constraint_groups);
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
