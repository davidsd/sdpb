#include "../parse_simple_command_line.hxx"

#include <El.hpp>

// void read_input_files(
//   const std::vector<boost::filesystem::path> &input_files,
//   El::BigFloat &objective_const, std::vector<El::BigFloat> &dual_objective_b,
//   std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices);

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
      std::vector<boost::filesystem::path> input_files;
      boost::filesystem::path output_dir;

      parse_simple_command_line("sdp2blocks", argc, argv, precision,
                                input_files, output_dir);
      mpf_set_default_prec(precision);
      El::gmp::SetPrecision(precision);
      El::mpfr::SetPrecision(precision);

      // El::BigFloat objective_const;
      // std::vector<El::BigFloat> dual_objective_b;
      // std::vector<Dual_Constraint_Group> dual_constraint_groups;
      // {
      //   std::vector<Polynomial_Vector_Matrix> polynomial_vector_matrices;
      //   read_input_files(input_files, objective_const, dual_objective_b,
      //                    polynomial_vector_matrices);
      //   for(auto &m : polynomial_vector_matrices)
      //     {
      //       dual_constraint_groups.emplace_back(m);
      //     }
      // }

      // boost::filesystem::create_directories(output_dir);
      // write_objectives(output_dir, objective_const, dual_objective_b);
      // write_bilinear_bases(output_dir, dual_constraint_groups);
      // write_blocks(output_dir, dual_constraint_groups);
      // write_primal_objective_c(output_dir, dual_constraint_groups);
      // write_free_var_matrix(output_dir, dual_objective_b.size(),
      //                       dual_constraint_groups);
    }
  catch(std::runtime_error &e)
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
