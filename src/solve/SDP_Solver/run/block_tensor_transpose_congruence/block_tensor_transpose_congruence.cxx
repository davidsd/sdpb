#include "../../../SDP_Solver.hxx"

// result[b] = Q[b]'^T A[b] Q[b]' for each block 0 <= b < Q.size()
// result[b], A[b] denote the b-th blocks of result, A, resp.

// Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product

// for each b, L.blocks[b], Q[b], Work[b], and result.blocks[b] must
// have the structure described above for `tensorTransposeCongruence'

void tensor_transpose_congruence(const Matrix &A, const Matrix &base,
                                 Matrix &Work, Matrix &result);

void block_tensor_transpose_congruence(
  const Block_Diagonal_Matrix &A, const std::vector<Matrix> &bilinear_bases,
  std::vector<Matrix> &Work, Block_Diagonal_Matrix &result)
{
  for(unsigned int b = 0; b < bilinear_bases.size(); b++)
    tensor_transpose_congruence(A.blocks[b], bilinear_bases[b], Work[b],
                                result.blocks[b]);
}

void block_tensor_transpose_congruence(
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result)
{
  for(unsigned int b = 0; b < bilinear_bases.size(); b++)
    {
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
           El::BigFloat(1), Y.blocks_elemental[b], bilinear_bases[b], 
           El::BigFloat(0), workspace[b]);
      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), bilinear_bases[b], workspace[b], 
           El::BigFloat(0), result.blocks_elemental[b]);
    }
}
