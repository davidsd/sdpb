#include "../../../Block_Diagonal_Matrix.hxx"

// A_X_inv = bilinear_base^T X^{-1} bilinear_base for each block

void compute_A_X_inv(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &A_X_inv)
{
  auto X_cholesky_block(X_cholesky.blocks.begin());
  auto A_X_inv_block(A_X_inv.blocks.begin());
  auto bilinear_bases_block(bilinear_bases.begin());

  for(auto &work : workspace)
    {
      // FIXME: This should be a constant, not computed over and over.
      // Set up the workspace to have copies of bilinear_bases[b]
      // along the diagonal
      for(int64_t row = 0; row < work.LocalHeight(); ++row)
        {
          size_t global_row(work.GlobalRow(row));
          size_t row_block(global_row / bilinear_bases_block->Height());

          for(int64_t column = 0; column < work.LocalWidth(); ++column)
            {
              size_t global_column(work.GlobalCol(column));
              size_t column_block(global_column
                                  / bilinear_bases_block->Width());
              work.SetLocal(
                row, column,
                row_block != column_block
                  ? El::BigFloat(0)
                  : (*bilinear_bases_block)(
                      global_row % bilinear_bases_block->Height(),
                      global_column % bilinear_bases_block->Width()));
            }
        }

      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), *X_cholesky_block, work);

      // We have to set this to zero because the values can be NaN.
      // Multiplying 0*NaN = NaN.
      // FIXME: GMP does not have NaN.  So remove?
      Zero(*A_X_inv_block);
      Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE,
           El::BigFloat(1), work, El::BigFloat(0), *A_X_inv_block);
      El::MakeSymmetric(El::UpperOrLower::LOWER, *A_X_inv_block);
      ++X_cholesky_block;
      ++A_X_inv_block;
      ++bilinear_bases_block;
    }
}
