#include "../../../Block_Diagonal_Matrix.hxx"

// Q_X_inv_Q = bilinear_base^T X^{-1} bilinear_base for each block

void compute_Q_X_inv_Q(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &Q_X_inv_Q)
{
  auto X_cholesky_block(X_cholesky.blocks.begin());
  auto Q_X_inv_Q_block(Q_X_inv_Q.blocks.begin());

  for(auto &block : bases_blocks)
    {
      auto temp_space(block);
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), *X_cholesky_block, temp_space);

      // We have to set this to zero because the values can be NaN.
      // Multiplying 0*NaN = NaN.
      Zero(*Q_X_inv_Q_block);
      Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE,
           El::BigFloat(1), temp_space, El::BigFloat(0), *Q_X_inv_Q_block);
      El::MakeSymmetric(El::UpperOrLower::LOWER, *Q_X_inv_Q_block);
      ++X_cholesky_block;
      ++Q_X_inv_Q_block;
    }
}
