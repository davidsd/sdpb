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
  std::array<std::vector<std::vector<std::vector<El::Matrix<El::BigFloat>>>>,
             2> &Q_Y_Q)
{
  Q_Y_Q[0].resize(bases_blocks.size() / 2);
  Q_Y_Q[1].resize(bases_blocks.size() / 2);

  for(size_t index(0); index < bases_blocks.size(); ++index)
    {
      auto &block(bases_blocks[index]);
      auto &Y_block(Y.blocks[index]);

      El::DistMatrix<El::BigFloat> Y_Q(block),
        Q_Y_Q_matrix_temp(block.Width(), block.Width(), block.Grid());
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
           Y_block, block, El::BigFloat(0), Y_Q);

      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), block, Y_Q, El::BigFloat(0), Q_Y_Q_matrix_temp);

      const size_t block_size(
        block_info.num_points.at(block_info.block_indices.at(index / 2))),
        dim(block_info.dimensions.at(block_info.block_indices.at(index / 2))),
        parity(index % 2), Q_index(index / 2);
      auto &Q_Y_Q_block(Q_Y_Q[parity][Q_index]);
      Q_Y_Q_block.resize(dim);

      El::MakeSymmetric(El::UpperOrLower::LOWER, Q_Y_Q_matrix_temp);

      El::DistMatrix<El::BigFloat,El::STAR,El::STAR> Y_star(block_size, block_size, block.Grid());
      for(size_t column_block = 0; column_block < dim; ++column_block)
        {
          Q_Y_Q_block[column_block].clear();
          Q_Y_Q_block[column_block].reserve(dim);
          const size_t column_offset(column_block * block_size);
          for(size_t row_block = 0; row_block < dim; ++row_block)
            {
              const size_t row_offset(row_block * block_size);
              El::DistMatrix<El::BigFloat> submatrix(
                El::LockedView(Q_Y_Q_matrix_temp, column_offset, row_offset,
                               block_size, block_size));

              Q_Y_Q_block[column_block].emplace_back(block_size, block_size);
              auto &Y_local(Q_Y_Q_block[column_block].back());
              
              El::Transpose(submatrix, Y_star);
              for(int64_t row(0); row<Y_star.LocalHeight(); ++row)
                for(int64_t column(0); column<Y_star.LocalWidth(); ++column)
                  {
                    Y_local(row,column)=Y_star.GetLocal(row,column);
                  }
            }
        }
    }
}
