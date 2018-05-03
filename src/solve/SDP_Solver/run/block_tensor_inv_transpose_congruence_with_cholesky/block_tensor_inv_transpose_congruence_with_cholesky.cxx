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
  for(size_t b = 0; b < bilinear_bases.size(); b++)
    {
      tensor_inv_transpose_congruence_with_cholesky(
        X_cholesky.blocks[b], bilinear_bases[b], workspace[b],
        result.blocks[b]);
    }
}

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result)
{
  for(size_t b = 0; b < bilinear_bases.size(); b++)
    {
      // Set up the workspace[b] to have copies of bilinear_bases[b]
      // along the diagonal
      size_t row_offset(workspace[b].GlobalRow(0)),
        column_offset(workspace[b].GlobalCol(0));

      for(int64_t row = 0; row < workspace[b].LocalHeight(); ++row)
        for(int64_t column = 0; column < workspace[b].LocalWidth(); ++column)
          {
            size_t m_row((row + row_offset) / bilinear_bases[b].Height()),
              m_column((column + column_offset) / bilinear_bases[b].Width());
            workspace[b].SetLocal(
              row, column,
              m_row != m_column
                ? El::BigFloat(0)
                : bilinear_bases[b].Get(
                    (row + row_offset) % bilinear_bases[b].Height(),
                    (column + column_offset) % bilinear_bases[b].Width()));
          }
      El::cholesky::SolveAfter(El::UpperOrLowerNS::LOWER,
                               El::Orientation::NORMAL,
                               X_cholesky.blocks_elemental[b], workspace[b]);

      Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE,
           El::BigFloat(1), workspace[b], El::BigFloat(0),
           result.blocks_elemental[b]);
    }
}
