#include "../../SDP_Solver.hxx"

// result[b] = Q[b]'^T A[b] Q[b]' for each block 0 <= b < Q.size()
// result[b], A[b] denote the b-th blocks of result, A, resp.

// Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product

// for each b, L.blocks[b], Q[b], Work[b], and result.blocks[b] must
// have the structure described above for `tensorTransposeCongruence'

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
      for(int64_t row = 0; row < workspace[b].LocalHeight(); ++row)
        {
          size_t global_row(workspace[b].GlobalRow(row));

          for(int64_t column = 0; column < workspace[b].LocalWidth(); ++column)
            {
              size_t global_column(workspace[b].GlobalCol(column));

              size_t row_block(global_row / bilinear_bases[b].Height()),
                column_block(global_column / bilinear_bases[b].Width());
              workspace[b].SetLocal(
                row, column,
                row_block != column_block
                  ? El::BigFloat(0)
                  : bilinear_bases[b].Get(
                      global_row % bilinear_bases[b].Height(),
                      global_column % bilinear_bases[b].Width()));
            }
        }
      auto temp_space(workspace[b]);
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
           Y.blocks_elemental[b], workspace[b], El::BigFloat(0), temp_space);
      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), workspace[b], temp_space, El::BigFloat(0),
           result.blocks_elemental[b]);
    }
}
