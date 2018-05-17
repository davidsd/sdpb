#include "../Block_Diagonal_Matrix.hxx"
#include "../Block_Matrix.hxx"

// B := L^{-1} B, where L is lower-triangular
void lower_triangular_solve(const Block_Diagonal_Matrix &L, Matrix &B)
{
}

// v := L^{-1} v, where L is lower-triangular
void lower_triangular_solve(const Block_Diagonal_Matrix &L, Vector &v)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      lower_triangular_solve(L.blocks[b], &v[L.block_start_indices[b]], 1,
                             v.size());
    }
}

