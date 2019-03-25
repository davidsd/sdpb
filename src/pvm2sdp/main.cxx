#include "../sdp_convert.hxx"

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<boost::filesystem::path> &input_files,
                        boost::filesystem::path &output_dir);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files, const int &rank,
  const int &num_procs, El::BigFloat &objective_const,
  std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Polynomial_Vector_Matrix> &local_polynomial_vector_matrices,
  std::vector<size_t> &indices);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  const int rank(El::mpi::Rank()),
    num_procs(El::mpi::Size(El::mpi::COMM_WORLD));

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
      {
        std::vector<Polynomial_Vector_Matrix> local_polynomial_vector_matrices;
        read_input_files(input_files, rank, num_procs, objective_const,
                         dual_objective_b, local_polynomial_vector_matrices,
                         indices);
        for(auto &pvm : local_polynomial_vector_matrices)
          {
            dual_constraint_groups.emplace_back(pvm);
          }
      }
      write_sdpb_input_files(output_dir, rank, num_procs, indices,
                             objective_const, dual_objective_b,
                             dual_constraint_groups);
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
