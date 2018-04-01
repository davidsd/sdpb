#include "../Block_Diagonal_Matrix.hxx"

// Tr(A B), where A and B are symmetric
Real frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                 const Block_Diagonal_Matrix &B)
{
  Real result = 0;
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    {
      Real f = frobenius_product_symmetric(A.blocks[b], B.blocks[b]);
      // this pragma means that other threads should be stopped while
      // this operation is performed (this prevents result from being
      // modified by multiple threads simultaneously)
      {
        result += f;
      }
    }
  return result;
}
