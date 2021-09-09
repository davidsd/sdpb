#include "Format.hxx"
#include "../sdp_read.hxx"
#include "../sdp_convert.hxx"
#include "../sdp_solve.hxx"
#include "../read_vector.hxx"

#include <boost/filesystem.hpp>

void handle_arguments(const int &argc, char **argv, El::BigFloat &threshold,
                      Format &format, boost::filesystem::path &input_path,
                      boost::filesystem::path &solution_path,
                      boost::filesystem::path &output_path);

std::vector<El::Matrix<El::BigFloat>>
read_x(const boost::filesystem::path &solution_path,
       const std::vector<Polynomial_Vector_Matrix> &matrices);

std::vector<El::Matrix<El::BigFloat>>
read_x(const boost::filesystem::path &solution_path,
       const std::vector<Positive_Matrix_With_Prefactor> &matrices);

El::Matrix<El::BigFloat>
read_y(const boost::filesystem::path &solution_path, const size_t &y_height);

void write_spectrum(const boost::filesystem::path &output_path,
                    const std::vector<std::vector<El::BigFloat>> &zeros);

std::vector<std::vector<El::BigFloat>> compute_spectrum_pmp(
  const std::vector<El::BigFloat> &normalization,
  const El::Matrix<El::BigFloat> &y,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices,
  const El::BigFloat &threshold);

std::vector<std::vector<El::BigFloat>>
compute_spectrum_pvm(const El::Matrix<El::BigFloat> &y,
                     const std::vector<Polynomial_Vector_Matrix> &matrices,
                     const std::vector<El::Matrix<El::BigFloat>> &x,
                     const El::BigFloat &threshold);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      El::BigFloat threshold;
      Format format;
      boost::filesystem::path input_path, solution_dir, output_path;
      handle_arguments(argc, argv, threshold, format, input_path, solution_dir,
                       output_path);

      switch(format)
        {
          case Format::Polynomial_Vector_Matrix: {
            std::vector<El::BigFloat> objectives;
            std::vector<Polynomial_Vector_Matrix> matrices;
            size_t num_blocks(0);
            read_pvm_input({input_path}, objectives, matrices, num_blocks);
            El::Matrix<El::BigFloat> y(
              read_y(solution_dir, objectives.size() - 1));
            std::vector<El::Matrix<El::BigFloat>> x(
              read_x(solution_dir, matrices));
            const std::vector<std::vector<El::BigFloat>> zeros(
              compute_spectrum_pvm(y, matrices, x, threshold));
            write_spectrum(output_path, zeros);
          }
          break;
          case Format::Positive_Matrix_with_Prefactor: {
            std::vector<El::BigFloat> objectives, normalization;
            std::vector<Positive_Matrix_With_Prefactor> matrices;
            read_input(input_path, objectives, normalization, matrices);
            El::Matrix<El::BigFloat> y(
              read_y(solution_dir / "y.txt", objectives.size() - 1));
            std::vector<El::Matrix<El::BigFloat>> x(
              read_x(solution_dir, matrices));
            const std::vector<std::vector<El::BigFloat>> zeros(
              compute_spectrum_pmp(normalization, y, matrices, threshold));
            write_spectrum(output_path, zeros);
          }
          break;
        default: throw std::runtime_error("INTERNAL ERROR");
        }
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
