#include "../../SDP_Solver.hxx"

// result = bilinear_base^T X^{-1} bilinear_base for each block

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
      for(int64_t row = 0; row < workspace[b].LocalHeight(); ++row)
        {
          size_t global_row(workspace[b].GlobalRow(row));
          size_t row_block(global_row / bilinear_bases[b].Height());

          for(int64_t column = 0; column < workspace[b].LocalWidth(); ++column)
            {
              size_t global_column(workspace[b].GlobalCol(column));
              size_t column_block(global_column / bilinear_bases[b].Width());

              workspace[b].SetLocal(
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
               El::BigFloat(1), X_cholesky.blocks_elemental[b], workspace[b]);

      // We have to set this to zero because the values can be NaN.
      // Multiplying 0*NaN = NaN.
      Zero(result.blocks_elemental[b]);
      Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE,
           El::BigFloat(1), workspace[b], El::BigFloat(0),
           result.blocks_elemental[b]);
      El::MakeSymmetric(El::UpperOrLower::LOWER, result.blocks_elemental[b]);
    }
}
