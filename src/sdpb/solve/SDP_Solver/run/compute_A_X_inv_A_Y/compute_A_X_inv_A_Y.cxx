#include "../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../Timers.hxx"

void compute_A_X_inv(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &A_X_inv);

void compute_A_Y(const Block_Diagonal_Matrix &Y,
                 const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
                 std::vector<El::DistMatrix<El::BigFloat>> &workspace,
                 Block_Diagonal_Matrix &A_Y);

void compute_A_X_inv_A_Y(
  const Block_Diagonal_Matrix &X_cholesky, const Block_Diagonal_Matrix &Y,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &A_X_inv, Block_Diagonal_Matrix &A_Y, Timers &timers)
{
  auto &congruence_timer(timers.add_and_start("run.bilinear_pairings"));
  compute_A_X_inv(X_cholesky, bilinear_bases, workspace, A_X_inv);

  compute_A_Y(Y, bilinear_bases, workspace, A_Y);
  congruence_timer.stop();
}
