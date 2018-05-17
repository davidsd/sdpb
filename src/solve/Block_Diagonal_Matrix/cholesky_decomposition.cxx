#include "../Block_Diagonal_Matrix.hxx"

// Compute L (lower triangular) such that A = L L^T
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      // FIXME: Use pivoting?
      //
      // Note: This only sets the lower triangular part and diagonal
      // of L.  The upper triangular part of L is still equal to A.
      L.blocks_elemental[b] = A.blocks_elemental[b];
      Cholesky(El::UpperOrLowerNS::LOWER, L.blocks_elemental[b]);
    }
}
