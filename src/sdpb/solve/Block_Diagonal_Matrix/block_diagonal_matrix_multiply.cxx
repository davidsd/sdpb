#include "../Block_Diagonal_Matrix.hxx"

// C := A*B
void block_diagonal_matrix_multiply(const Block_Diagonal_Matrix &A,
                                    const Block_Diagonal_Matrix &B,
                                    Block_Diagonal_Matrix &C)
{
  block_diagonal_matrix_scale_multiply_add(El::BigFloat(1), A, B,
                                           El::BigFloat(0), C);
}
