#include "../Block_Diagonal_Matrix.hxx"

// C := alpha*A*B + beta*C
void block_diagonal_matrix_scale_multiply_add(const El::BigFloat &alpha,
                                              const Block_Diagonal_Matrix &A,
                                              const Block_Diagonal_Matrix &B,
                                              const El::BigFloat &beta,
                                              Block_Diagonal_Matrix &C)
{
  for(size_t block = 0; block < A.blocks_elemental.size(); ++block)
    {
      El::Gemm(El::OrientationNS::NORMAL, El::OrientationNS::NORMAL, alpha,
               A.blocks_elemental[block], B.blocks_elemental[block], beta,
               C.blocks_elemental[block]);
    }
}
