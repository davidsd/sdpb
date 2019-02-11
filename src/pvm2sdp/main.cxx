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

  if(El::mpi::Size(El::mpi::COMM_WORLD) != 1)
    {
      if(El::mpi::Rank() == 0)
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
        // FIXME: Make this work with multiple MPI processes.
        for(auto &m : polynomial_vector_matrices)
          {
            dual_constraint_groups.emplace_back(El::mpi::Rank(), m);
          }
      }

      boost::filesystem::create_directories(output_dir);
      write_objectives(output_dir, objective_const, dual_objective_b);
      write_bilinear_bases(output_dir, dual_constraint_groups);
      write_blocks(output_dir, dual_constraint_groups);
      write_primal_objective_c(output_dir, dual_constraint_groups);
      write_free_var_matrix(output_dir, dual_objective_b.size(),
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
