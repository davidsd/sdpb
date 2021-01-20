#include "../../../Block_Diagonal_Matrix.hxx"
#include "../../../../Timers.hxx"

void compute_Q_X_inv_Q(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &Q_X_inv_Q);

void compute_Q_Y_Q(
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &Q_Y_Q);

void compute_bilinear_pairings(
  const Block_Diagonal_Matrix &X_cholesky, const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &Q_X_inv_Q, Block_Diagonal_Matrix &Q_Y_Q,
  Timers &timers)
{
  auto &congruence_timer(timers.add_and_start("run.bilinear_pairings"));
  compute_Q_X_inv_Q(X_cholesky, bases_blocks, Q_X_inv_Q);

  compute_Q_Y_Q(Y, bases_blocks, Q_Y_Q);
  congruence_timer.stop();
}
