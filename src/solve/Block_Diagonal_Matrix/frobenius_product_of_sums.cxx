#include "../Block_Diagonal_Matrix.hxx"

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
El::BigFloat
frobenius_product_of_sums_elemental(const Block_Diagonal_Matrix &X,
                                    const Block_Diagonal_Matrix &dX,
                                    const Block_Diagonal_Matrix &Y,
                                    const Block_Diagonal_Matrix &dY)
{
  // FIXME: This can be sped up by not have intermediate results.  It
  // may require looking into the implementation of Dotu.
  El::BigFloat result = 0;
  for(size_t b = 0; b < X.blocks_elemental.size(); b++)
    {
      El::DistMatrix<El::BigFloat> X_dX(X.blocks_elemental[b]);
      X_dX += dX.blocks_elemental[b];
      El::DistMatrix<El::BigFloat> Y_dY(Y.blocks_elemental[b]);
      Y_dY += dY.blocks_elemental[b];
      result += Dotu(X_dX, Y_dY);
    }
  return result;
}
