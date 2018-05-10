#include "../Block_Diagonal_Matrix.hxx"

// C := alpha*A*B + beta*C
void block_diagonal_matrix_scale_multiply_add(const Real &alpha,
                                              Block_Diagonal_Matrix &A,
                                              Block_Diagonal_Matrix &B,
                                              const Real &beta,
                                              Block_Diagonal_Matrix &C)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      matrix_scale_multiply_add(alpha, A.blocks[b], B.blocks[b], beta,
                                C.blocks[b]);
    }
}

void block_diagonal_matrix_scale_multiply_add(const El::BigFloat &alpha,
                                              const Block_Diagonal_Matrix &A,
                                              const Block_Diagonal_Matrix &B,
                                              const El::BigFloat &beta,
                                              Block_Diagonal_Matrix &C)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      El::Gemm(El::OrientationNS::NORMAL, El::OrientationNS::NORMAL, alpha,
               A.blocks_elemental[b], B.blocks_elemental[b], beta,
               C.blocks_elemental[b]);
    }
}
