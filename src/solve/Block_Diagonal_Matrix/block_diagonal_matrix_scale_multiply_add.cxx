#include "../Block_Diagonal_Matrix.hxx"

// C := alpha*A*B + beta*C
void block_diagonal_matrix_scale_multiply_add(Real alpha,
                                              Block_Diagonal_Matrix &A,
                                              Block_Diagonal_Matrix &B,
                                              Real beta,
                                              Block_Diagonal_Matrix &C)
{
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    {
      matrix_scale_multiply_add(alpha, A.blocks[b], B.blocks[b], beta,
                                C.blocks[b]);
    }
}
