#include "sdp_solve/sdp_solve.hxx"

#include <filesystem>

namespace fs = std::filesystem;

// TODO: Have this be part of sdp_solve.hxx
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);
void compute_A_X_inv(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_X_inv);

void compute_A_Y(
  const Block_Info &block_info, const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_Y);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const El::Grid &block_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q,
  Timers &timers);

void setup_solver(const Block_Info &block_info, const El::Grid &grid,
                  const SDP &sdp, const fs::path &solution_dir,
                  Block_Diagonal_Matrix &schur_complement_cholesky,
                  Block_Matrix &schur_off_diagonal,
                  El::DistMatrix<El::BigFloat> &Q)
{
  if(fs::exists(solution_dir / "Q_cholesky.txt"))
    {
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          read_text_block(schur_complement_cholesky.blocks.at(block),
                          solution_dir, "schur_complement_cholesky_",
                          block_index);
          read_text_block(schur_off_diagonal.blocks.at(block),
                          solution_dir, "schur_off_diagonal_", block_index);
        }
      read_text_block(Q, solution_dir / "Q_cholesky.txt");
    }
  else
    {
      Block_Diagonal_Matrix X(block_info.psd_matrix_block_sizes(),
                              block_info.block_indices,
                              block_info.num_points.size(), grid),
        Y(X);
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          for(size_t psd_block(0); psd_block < 2; ++psd_block)
            {
              // Constant constraints have empty odd parity blocks, so we do
              // not need to load them.
              if(X.blocks.at(2 * block + psd_block).Height() != 0)
                {
                  const size_t psd_index(2 * block_index + psd_block);
                  read_text_block(X.blocks.at(2 * block + psd_block),
                                  solution_dir, "X_matrix_", psd_index);
                  read_text_block(Y.blocks.at(2 * block + psd_block),
                                  solution_dir, "Y_matrix_", psd_index);
                }
            }
        }

      std::array<
        std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
        A_X_inv, A_Y;

      Block_Diagonal_Matrix X_cholesky(X);
      cholesky_decomposition(X, X_cholesky);
      compute_A_X_inv(block_info, X_cholesky, sdp.bases_blocks, A_X_inv);
      compute_A_Y(block_info, Y, sdp.bases_blocks, A_Y);

      Timers timers;
      initialize_schur_complement_solver(block_info, sdp, A_X_inv, A_Y, grid,
                                         schur_complement_cholesky,
                                         schur_off_diagonal, Q, timers);
    }
}
