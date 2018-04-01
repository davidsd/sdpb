#include "../Block_Diagonal_Matrix.hxx"

// Compute L (lower triangular) such that A + U U^T = L L^T, where U
// is a low-rank update matrix.
void cholesky_decomposition_stabilized(
  Block_Diagonal_Matrix &A, Block_Diagonal_Matrix &L,
  vector<vector<Integer>> &stabilize_indices,
  vector<vector<Real>> &stabilize_lambdas, const double stabilize_threshold)
{
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    {
      cholesky_decomposition_stabilized(
        A.blocks[b], L.blocks[b], stabilize_indices[b], stabilize_lambdas[b],
        stabilize_threshold);
    }
}
