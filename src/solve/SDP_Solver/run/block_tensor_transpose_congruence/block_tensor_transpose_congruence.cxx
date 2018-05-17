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

// FIXME: Why am I using the DistMatrix version of bilinear_bases
// instead of the local version?  With the local version, there would
// be no communication required.
void block_tensor_transpose_congruence(
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result)
{
  for(size_t b = 0; b < bilinear_bases.size(); b++)
    {
      // FIXME: This should be a constant, not computed over and over.
      // Set up the workspace[b] to have copies of bilinear_bases[b]
      // along the diagonal
      size_t row_offset(workspace[b].GlobalRow(0)),
        column_offset(workspace[b].GlobalCol(0));

      for(int64_t row = 0; row < workspace[b].LocalHeight(); ++row)
        for(int64_t column = 0; column < workspace[b].LocalWidth(); ++column)
          {
            size_t m_row((row + row_offset) / bilinear_bases[b].Height()),
              m_column((column + column_offset) / bilinear_bases[b].Width());
            workspace[b].SetLocal(
              row, column,
              m_row != m_column
                ? El::BigFloat(0)
                : bilinear_bases[b].Get(
                    (row + row_offset) % bilinear_bases[b].Height(),
                    (column + column_offset) % bilinear_bases[b].Width()));
          }
      
      auto temp_space(workspace[b]);
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
           El::BigFloat(1), Y.blocks_elemental[b], workspace[b], 
           El::BigFloat(0), temp_space);
      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), workspace[b], temp_space, 
           El::BigFloat(0), result.blocks_elemental[b]);
    }
}
