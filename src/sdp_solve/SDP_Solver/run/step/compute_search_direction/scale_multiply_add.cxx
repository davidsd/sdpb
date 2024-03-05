#include "sdp_solve/Block_Diagonal_Matrix.hxx"

// C := alpha*A*B + beta*C
void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C)
{
  for(size_t block = 0; block < A.blocks.size(); ++block)
    {
      El::Gemm(El::OrientationNS::NORMAL, El::OrientationNS::NORMAL, alpha,
               A.blocks[block], B.blocks[block], beta, C.blocks[block]);
    }
}
