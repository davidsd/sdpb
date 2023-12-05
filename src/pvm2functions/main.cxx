#include "sdpb_util/Boost_Float.hxx"
#include "sdp_read/sdp_read.hxx"

namespace fs = std::filesystem;

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<fs::path> &input_files,
                        fs::path &output_dir);

void write_functions(
  const fs::path &output_path,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      int precision;
      std::vector<fs::path> input_files;
      fs::path output_path;

      parse_command_line(argc, argv, precision, input_files, output_path);
      El::gmp::SetPrecision(precision);
      Boost_Float::default_precision(precision * log(2) / log(10));

      std::vector<El::BigFloat> dual_objective_b;
      std::vector<Polynomial_Vector_Matrix> polynomial_vector_matrices;
      size_t num_blocks(0);
      read_pvm_input(input_files, dual_objective_b, polynomial_vector_matrices,
                     num_blocks);

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
