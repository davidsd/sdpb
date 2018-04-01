#include "../Block_Diagonal_Matrix.hxx"

// Tr(A B), where A and B are symmetric
Real frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                 const Block_Diagonal_Matrix &B)
{
  Real result = 0;
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      result += frobenius_product_symmetric(A.blocks[b], B.blocks[b]);
    }
  return result;
}
