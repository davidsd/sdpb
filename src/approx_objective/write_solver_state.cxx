#include "sdp_solve/sdp_solve.hxx"
#include "sdpb_util/write_distmatrix.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void write_solver_state(const std::vector<size_t> &block_indices,
                        const fs::path &solution_dir,
                        const Block_Diagonal_Matrix &schur_complement_cholesky,
                        const Block_Matrix &schur_off_diagonal,
                        const El::DistMatrix<El::BigFloat> &Q)
{
  for(size_t block = 0; block != schur_complement_cholesky.blocks.size();
      ++block)
    {
      size_t block_index(block_indices.at(block));
      write_distmatrix(schur_complement_cholesky.blocks.at(block),
                       solution_dir
                         / ("schur_complement_cholesky_"
                            + std::to_string(block_index) + ".txt"));
      write_distmatrix(
        schur_off_diagonal.blocks.at(block),
        solution_dir
          / ("schur_off_diagonal_" + std::to_string(block_index) + ".txt"));
    }
  write_distmatrix(Q, solution_dir / "Q_cholesky.txt");
}
