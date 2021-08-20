#include "../sdp_read.hxx"
#include "../sdp_solve.hxx"
#include "../read_vector.hxx"

#include <boost/filesystem.hpp>

void handle_arguments(const int &argc, char **argv,
                      El::BigFloat &threshold,
                      boost::filesystem::path &input_path,
                      boost::filesystem::path &solution_path,
                      boost::filesystem::path &output_path);

void write_spectrum(const boost::filesystem::path &output_path,
                    const std::vector<std::vector<El::BigFloat>> &zeros);

std::vector<std::vector<El::BigFloat>>
compute_spectrum(const std::vector<El::BigFloat> &normalization,
                 const El::DistMatrix<El::BigFloat> &y,
                 const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                 const El::BigFloat &threshold);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      El::BigFloat threshold;
      boost::filesystem::path input_path, solution_path, output_path;
      handle_arguments(argc, argv, threshold, input_path, solution_path, output_path);

      std::vector<El::BigFloat> objectives, normalization;
      std::vector<Positive_Matrix_With_Prefactor> matrices;
      read_input(input_path, objectives, normalization, matrices);
      El::DistMatrix<El::BigFloat> y(normalization.size() - 1, 1);
      read_text_block(y, solution_path);
      const std::vector<std::vector<El::BigFloat>> zeros(compute_spectrum(
        normalization, y, matrices, threshold));
      write_spectrum(output_path, zeros);
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
