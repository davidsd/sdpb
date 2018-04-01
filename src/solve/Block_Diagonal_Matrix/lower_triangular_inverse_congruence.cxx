#include "../Block_Diagonal_Matrix.hxx"

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(Block_Diagonal_Matrix &A,
                                         Block_Diagonal_Matrix &L)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      lower_triangular_inverse_congruence(A.blocks[b], L.blocks[b]);
    }
}
