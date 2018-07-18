#include "Dual_Constraint_Group.hxx"

#include <boost/filesystem.hpp>
#include <vector>

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<boost::filesystem::path> &input_files,
                        boost::filesystem::path &output_dir);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  std::vector<El::BigFloat> &objective,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices);

void write_objective(const boost::filesystem::path &output_dir,
                     const std::vector<El::BigFloat> &objective);

void write_bilinear_bases(
  const boost::filesystem::path &output_dir,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  int result(0);
  try
    {
      int precision;
      std::vector<boost::filesystem::path> input_files;
      boost::filesystem::path output_dir;

      parse_command_line(argc, argv, precision, input_files, output_dir);
      mpf_set_default_prec(precision);
      El::gmp::SetPrecision(precision);
      El::mpfr::SetPrecision(precision);

      std::vector<El::BigFloat> objective;
      std::vector<Dual_Constraint_Group> dual_constraint_groups;
      {
        std::vector<Polynomial_Vector_Matrix> polynomial_vector_matrices;
        read_input_files(input_files, objective, polynomial_vector_matrices);
        for(auto &m : polynomial_vector_matrices)
          {
            dual_constraint_groups.emplace_back(m);
          }
      }

      boost::filesystem::create_directories(output_dir);
      write_objective(output_dir, objective);
      write_bilinear_bases(output_dir, dual_constraint_groups);
      // write_blocks(output_dir, dual_constraint_groups);

      // write_primal_objective(output_dir, dual_constraint_groups);
      // write_free_var_matrix(output_dir, dual_constraint_groups);
    }
  catch(std::runtime_error &e)
    {
      std::cerr << "Error: " << e.what() << "\n" << std::flush;
      result = 1;
    }
  catch(...)
    {
      std::cerr << "Unknown Error\n" << std::flush;
      result = 1;
    }
  return result;
}
