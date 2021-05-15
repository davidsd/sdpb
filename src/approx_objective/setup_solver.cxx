#include "../sdp_solve.hxx"

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
                  const SDP &sdp,
                  const Block_Diagonal_Matrix &X,
                  const Block_Diagonal_Matrix &Y,
                  Block_Diagonal_Matrix &schur_complement_cholesky,
                  Block_Matrix &schur_off_diagonal,
                  El::DistMatrix<El::BigFloat> &Q)
{
  std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    A_X_inv, A_Y;

  Block_Diagonal_Matrix X_cholesky(X);
  cholesky_decomposition(X, X_cholesky);
  compute_A_X_inv(block_info, X_cholesky, sdp.bases_blocks, A_X_inv);
  compute_A_Y(block_info, Y, sdp.bases_blocks, A_Y);

  Timers timers(false);
  initialize_schur_complement_solver(block_info, sdp, A_X_inv, A_Y, grid,
                                     schur_complement_cholesky,
                                     schur_off_diagonal, Q, timers);
}
