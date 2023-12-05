#include "sdp_solve/Block_Diagonal_Matrix.hxx"

// Compute L (lower triangular) such that A = L L^T
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      // FIXME: Use pivoting?
      L.blocks[b] = A.blocks[b];
      Cholesky(El::UpperOrLowerNS::LOWER, L.blocks[b]);
    }
}
