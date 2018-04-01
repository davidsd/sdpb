#include "../Block_Diagonal_Matrix.hxx"

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void block_matrix_solve_with_cholesky(Block_Diagonal_Matrix &ACholesky,
                                      Block_Diagonal_Matrix &X)
{
  for(size_t b = 0; b < X.blocks.size(); b++)
    {
      matrix_solve_with_cholesky(ACholesky.blocks[b], X.blocks[b]);
    }
}
