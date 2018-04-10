#include "../Block_Diagonal_Matrix.hxx"

// B := L^{-1} B, where L is lower-triangular
void block_matrix_lower_triangular_solve(const Block_Diagonal_Matrix &L,
                                         Matrix &B)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      lower_triangular_solve(L.blocks[b], &B.elt(L.block_start_indices[b], 0),
                             B.cols, B.rows);
    }
}

// v := L^{-1} v, where L is lower-triangular
void block_matrix_lower_triangular_solve(const Block_Diagonal_Matrix &L,
                                         Vector &v)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      lower_triangular_solve(L.blocks[b], &v[L.block_start_indices[b]], 1,
                             v.size());
    }
}
