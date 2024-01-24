#include "Zeros.hxx"
#include "pmp2sdp/write_sdp.hxx"
#include "pmp_read/pmp_read.hxx"
#include "sdp_solve/sdp_solve.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void handle_arguments(const int &argc, char **argv, El::BigFloat &threshold,
                      El::BigFloat &mesh_threshold, fs::path &input_path,
                      fs::path &solution_path, fs::path &output_path,
                      bool &need_lambda);

std::vector<El::Matrix<El::BigFloat>>
read_x(const fs::path &solution_path,
       const std::vector<Polynomial_Vector_Matrix> &matrices,
       const std::vector<size_t> &block_indices);

El::Matrix<El::BigFloat>
read_y(const fs::path &solution_path, const size_t &y_height);

void write_spectrum(const fs::path &output_path, const size_t &num_blocks,
                    const std::vector<Zeros> &zeros_blocks,
                    const std::vector<size_t> &block_indices);

std::vector<Zeros>
compute_spectrum(const std::vector<El::BigFloat> &normalization,
                 const El::Matrix<El::BigFloat> &y,
                 const std::vector<Polynomial_Vector_Matrix> &matrices,
                 const std::vector<El::Matrix<El::BigFloat>> &x,
                 const El::BigFloat &threshold,
                 const El::BigFloat &mesh_threshold, const bool &need_lambda);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      El::BigFloat threshold, mesh_threshold;
      fs::path input_path, solution_dir, output_path;
      bool need_lambda;
      handle_arguments(argc, argv, threshold, mesh_threshold, input_path,
                       solution_dir, output_path, need_lambda);

      bool debug = false; // TODO set verbosity from command line
      Timers timers(env, debug);
      const auto pmp = read_polynomial_matrix_program(input_path, timers);
      const size_t num_blocks = pmp.num_matrices;
      const auto &block_indices = pmp.matrix_index_local_to_global;
      El::Matrix<El::BigFloat> y(pmp.objective.size() - 1, 1);
      read_text_block(y, solution_dir / "y.txt");
      std::vector<El::Matrix<El::BigFloat>> x(
        read_x(solution_dir, pmp.matrices, block_indices));
      const std::vector<Zeros> zeros_blocks(
        compute_spectrum(pmp.normalization, y, pmp.matrices, x, threshold,
                         mesh_threshold, need_lambda));
      write_spectrum(output_path, num_blocks, zeros_blocks, block_indices);
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
