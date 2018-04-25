#include "../Block_Diagonal_Matrix.hxx"
#include "../Block_Matrix.hxx"

// v := L^{-T} v, where L is lower-triangular
void block_matrix_lower_triangular_transpose_solve(
  const Block_Diagonal_Matrix &L, Vector &v)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      lower_triangular_transpose_solve(
        L.blocks[b], &v[L.block_start_indices[b]], 1, v.size());
    }
}

// v := L^{-T} v, where L is the result of a previous cholesky
// factorization
void block_matrix_lower_triangular_transpose_solve(
  const Block_Diagonal_Matrix &L, Block_Matrix &v)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      El::cholesky::SolveAfter(El::UpperOrLowerNS::UPPER,
                               El::OrientationNS::NORMAL,
                               L.blocks_elemental[b], v.blocks[b]);
    }
}
