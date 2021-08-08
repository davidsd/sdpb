#include "../Boost_Float.hxx"
#include "../sdp_convert.hxx"

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<boost::filesystem::path> &input_files,
                        boost::filesystem::path &output_dir);

void read_input_files(
  const std::vector<boost::filesystem::path> &input_files,
  std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices,
  size_t &num_processed);

void write_functions(
  const boost::filesystem::path &output_path,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      int precision;
      std::vector<boost::filesystem::path> input_files;
      boost::filesystem::path output_path;

      parse_command_line(argc, argv, precision, input_files, output_path);
      El::gmp::SetPrecision(precision);
      Boost_Float::default_precision(precision * log(2) / log(10));
      
      std::vector<El::BigFloat> dual_objective_b;
      std::vector<Polynomial_Vector_Matrix> polynomial_vector_matrices;
      size_t num_blocks(0);
      read_input_files(input_files, dual_objective_b,
                       polynomial_vector_matrices, num_blocks);

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
      write_functions(output_path, dual_objective_b,
                      polynomial_vector_matrices);
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
