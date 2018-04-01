#include "../Block_Diagonal_Matrix.hxx"

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
Real frobenius_product_of_sums(const Block_Diagonal_Matrix &X,
                               const Block_Diagonal_Matrix &dX,
                               const Block_Diagonal_Matrix &Y,
                               const Block_Diagonal_Matrix &dY)
{
  Real result = 0;
  for(unsigned int b = 0; b < X.blocks.size(); b++)
    {
      Real f = frobenius_product_of_sums(X.blocks[b], dX.blocks[b],
                                         Y.blocks[b], dY.blocks[b]);
      {
        result += f;
      }
    }
  return result;
}
