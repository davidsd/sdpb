#include "../../../SDP_Solver.hxx"

// result = bilinear_base^T X^{-1} bilinear_base for each block

void tensor_inv_transpose_congruence_with_cholesky(
  const Matrix &X_cholesky_block, const Matrix &bilinear_base_block,
  Matrix &workspace_block, Matrix &result);

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<Matrix> &bilinear_bases, std::vector<Matrix> &workspace,
  Block_Diagonal_Matrix &result)
{
  for(unsigned int b = 0; b < bilinear_bases.size(); b++)
    tensor_inv_transpose_congruence_with_cholesky(
      X_cholesky.blocks[b], bilinear_bases[b], workspace[b], result.blocks[b]);
}

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result)
{
  workspace = bilinear_bases;
  for(size_t b = 0; b < bilinear_bases.size(); b++)
    {
      El::cholesky::SolveAfter(
        El::UpperOrLowerNS::LOWER, El::Orientation::NORMAL,
        X_cholesky.blocks_elemental[b], workspace[b]);
      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), bilinear_bases[b], workspace[b],
           El::BigFloat(0), result.blocks_elemental[b]);
    }
}
