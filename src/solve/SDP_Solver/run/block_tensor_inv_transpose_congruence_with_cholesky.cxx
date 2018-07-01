#include "../../SDP_Solver.hxx"

// bilinear_pairings_X_inv = bilinear_base^T X^{-1} bilinear_base for each block

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  const std::vector<size_t> &block_indices,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &bilinear_pairings_X_inv)
{
  auto work(workspace.begin());
  auto X_cholesky_block(X_cholesky.blocks.begin());
  auto bilinear_pairings_X_inv_block(bilinear_pairings_X_inv.blocks.begin());
  for(auto &block_index : block_indices)
    {
      for(size_t b = 2 * block_index; b < 2 * block_index + 2; b++)
        {
          // Set up the workspace to have copies of bilinear_bases[b]
          // along the diagonal
          for(int64_t row = 0; row < work->LocalHeight(); ++row)
            {
              size_t global_row(work->GlobalRow(row));
              size_t row_block(global_row / bilinear_bases[b].Height());

              for(int64_t column = 0; column < work->LocalWidth(); ++column)
                {
                  size_t global_column(work->GlobalCol(column));
                  size_t column_block(global_column
                                      / bilinear_bases[b].Width());

                  work->SetLocal(
                    row, column,
                    row_block != column_block
                      ? El::BigFloat(0)
                      : bilinear_bases[b].Get(
                          global_row % bilinear_bases[b].Height(),
                          global_column % bilinear_bases[b].Width()));
                }
            }

          El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
                   El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
                   El::BigFloat(1), *X_cholesky_block, *work);

          // We have to set this to zero because the values can be NaN.
          // Multiplying 0*NaN = NaN.
          Zero(*bilinear_pairings_X_inv_block);
          Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE,
               El::BigFloat(1), *work, El::BigFloat(0),
               *bilinear_pairings_X_inv_block);
          El::MakeSymmetric(El::UpperOrLower::LOWER,
                            *bilinear_pairings_X_inv_block);
          ++work;
          ++X_cholesky_block;
          ++bilinear_pairings_X_inv_block;
        }
    }
}
