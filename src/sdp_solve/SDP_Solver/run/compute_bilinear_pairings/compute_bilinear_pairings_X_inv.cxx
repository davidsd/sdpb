#include "../../../Block_Diagonal_Matrix.hxx"

// bilinear_pairings_X_inv = bilinear_base^T X^{-1} bilinear_base for each block

void compute_bilinear_pairings_X_inv(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &bilinear_pairings_X_inv)
{
  auto X_cholesky_block(X_cholesky.blocks.begin());
  auto bilinear_pairings_X_inv_block(bilinear_pairings_X_inv.blocks.begin());

  for(auto &block : bases_blocks)
    {
      auto temp_space(block);
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), *X_cholesky_block, temp_space);

      // We have to set this to zero because the values can be NaN.
      // Multiplying 0*NaN = NaN.
      Zero(*bilinear_pairings_X_inv_block);
      Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE,
           El::BigFloat(1), temp_space, El::BigFloat(0),
           *bilinear_pairings_X_inv_block);
      El::MakeSymmetric(El::UpperOrLower::LOWER,
                        *bilinear_pairings_X_inv_block);
      ++X_cholesky_block;
      ++bilinear_pairings_X_inv_block;
    }
}
