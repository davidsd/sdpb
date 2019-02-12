#include "../sdp_convert.hxx"

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<boost::filesystem::path> &input_files,
                        boost::filesystem::path &output_dir);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objective_b,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  const int rank(El::mpi::Rank()),
    num_procs(El::mpi::Size(El::mpi::COMM_WORLD));
  if(num_procs != 1)
    {
      if(rank == 0)
        {
          std::cerr << "pvm2sdp can only be run with a single MPI task, but "
                       "was invoked with "
                    << El::mpi::Size(El::mpi::COMM_WORLD) << " tasks.\n"
                    << std::flush;
        }
      El::Finalize();
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }

  try
    {
      int precision;
      std::vector<boost::filesystem::path> input_files;
      boost::filesystem::path output_dir;

      parse_command_line(argc, argv, precision, input_files, output_dir);
      El::gmp::SetPrecision(precision);

      El::BigFloat objective_const;
      std::vector<El::BigFloat> dual_objective_b;
      std::vector<Dual_Constraint_Group> dual_constraint_groups;
      {
        std::vector<Polynomial_Vector_Matrix> polynomial_vector_matrices;
        read_input_files(input_files, objective_const, dual_objective_b,
                         polynomial_vector_matrices);
        for(auto &m : polynomial_vector_matrices)
          {
            dual_constraint_groups.emplace_back(m);
          }
      }
      std::vector<size_t> indices;
      for(size_t index = 0; index < dual_constraint_groups.size(); ++index)
        {
          indices.push_back(index);
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
