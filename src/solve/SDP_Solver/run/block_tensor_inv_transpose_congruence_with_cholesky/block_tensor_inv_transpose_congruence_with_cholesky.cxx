#include "../../../SDP_Solver.hxx"

// Result_b = Q[b]'^T A_b^{-1} Q[b]' for each block 0 <= b < Q.size()
// - Result_b, A_b denote the b-th blocks of Result, A, resp.
// - Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product
// - for each b, L.blocks[b], Q[b], Work[b], and Result.blocks[b] must
//   have the structure described above for
//   `tensorInvTransposeCongruenceWithCholesky'
//

void tensor_inv_transpose_congruence_with_cholesky(const Matrix &L,
                                                   const Matrix &Q,
                                                   Matrix &Work,
                                                   Matrix &Result);

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &L, const std::vector<Matrix> &Q,
  std::vector<Matrix> &Work, Block_Diagonal_Matrix &Result)
{
  for(unsigned int b = 0; b < Q.size(); b++)
    tensor_inv_transpose_congruence_with_cholesky(L.blocks[b], Q[b], Work[b],
                                                  Result.blocks[b]);
}
