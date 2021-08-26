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

void write_spectrum(const boost::filesystem::path &output_path,
                    const std::vector<std::vector<El::BigFloat>> &zeros);

std::vector<std::vector<El::BigFloat>>
compute_spectrum(const std::vector<El::BigFloat> &normalization,
                 const El::Matrix<El::BigFloat> &y,
                 const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                 const El::BigFloat &threshold);

std::vector<std::vector<El::BigFloat>>
compute_spectrum_pvm(const El::Matrix<El::BigFloat> &y,
                     const std::vector<Polynomial_Vector_Matrix> &matrices,
                     const El::BigFloat &threshold);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      El::BigFloat threshold;
      Format format;
      boost::filesystem::path input_path, solution_path, output_path;
      handle_arguments(argc, argv, threshold, format, input_path,
                       solution_path, output_path);

      switch(format)
        {
          case Format::Polynomial_Vector_Matrix: {
            std::vector<El::BigFloat> objectives;
            std::vector<Polynomial_Vector_Matrix> matrices;
            size_t num_blocks(0);
            read_pvm_input({input_path}, objectives, matrices, num_blocks);
            El::DistMatrix<El::BigFloat> y_dist(objectives.size() - 1, 1);
            read_text_block(y_dist, solution_path);
            El::DistMatrix<El::BigFloat, El::STAR, El::STAR> y_star(y_dist);
            const std::vector<std::vector<El::BigFloat>> zeros(
              compute_spectrum_pvm(y_star.LockedMatrix(), matrices,
                                   threshold));
            write_spectrum(output_path, zeros);
          }
          break;
          case Format::Positive_Matrix_with_Prefactor: {
            std::vector<El::BigFloat> objectives, normalization;
            std::vector<Positive_Matrix_With_Prefactor> matrices;
            read_input(input_path, objectives, normalization, matrices);
            El::DistMatrix<El::BigFloat> y_dist(normalization.size() - 1, 1);
            read_text_block(y_dist, solution_path);
            El::DistMatrix<El::BigFloat, El::STAR, El::STAR> y_star(y_dist);
            const std::vector<std::vector<El::BigFloat>> zeros(
              compute_spectrum(normalization, y_star.LockedMatrix(), matrices,
                               threshold));
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
