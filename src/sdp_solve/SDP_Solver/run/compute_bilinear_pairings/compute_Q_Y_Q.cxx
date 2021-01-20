#include "../../../Block_Diagonal_Matrix.hxx"
#include "../../../Block_Info.hxx"

// Q_Y_Q[b] = Q[b]'^T A[b] Q[b]' for each block 0 <= b < Q.size()
// Q_Y_Q[b], A[b] denote the b-th blocks of Q_Y_Q,
// A, resp.

// Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product

// for each b, L.blocks[b], Q[b], Work[b], and Q_Y_Q.blocks[b]
// must have the structure described above for `tensorTransposeCongruence'

void compute_Q_Y_Q(
  const Block_Info &block_info, const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &Q_Y_Q)
{
  for(size_t index(0); index < bases_blocks.size(); ++index)
    {
      auto &block(bases_blocks[index]);
      auto &Q_Y_Q_block(Q_Y_Q.blocks[index]);
      auto &Y_block(Y.blocks[index]);

      auto Y_Q(block), Q_Y_Q_temp(Q_Y_Q_block);
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
           Y_block, block, El::BigFloat(0), Y_Q);

      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), block, Y_Q, El::BigFloat(0), Q_Y_Q_temp);

      El::MakeSymmetric(El::UpperOrLower::LOWER, Q_Y_Q_temp);

      const size_t block_size(
        block_info.num_points.at(block_info.block_indices.at(index / 2)));
      El::DistMatrix<El::BigFloat> transpose(block_size, block_size,
                                             Q_Y_Q_temp.Grid());
      for(size_t column_block = 0; column_block < block_info.dimensions.at(
                                     block_info.block_indices.at(index / 2));
          ++column_block)
        {
          const size_t column_offset(column_block * block_size);
          for(size_t row_block = 0; row_block < block_info.dimensions.at(
                                      block_info.block_indices.at(index / 2));
              ++row_block)
            {
              const size_t row_offset(row_block * block_size);

              El::DistMatrix<El::BigFloat> submatrix(
                El::View(Q_Y_Q_temp, column_offset, row_offset, block_size,
                         block_size));
              El::DistMatrix<El::BigFloat> Q_Y_Q_submatrix(
                El::View(Q_Y_Q_block, column_offset, row_offset, block_size,
                         block_size));
              El::Transpose(submatrix, transpose);
              El::Copy(transpose, Q_Y_Q_submatrix);
            }
        }
    }
}
