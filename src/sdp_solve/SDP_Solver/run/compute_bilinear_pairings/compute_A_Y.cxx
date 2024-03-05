#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/Block_Info.hxx"

// A_Y[b] = Q[b]'^T A[b] Q[b]' for each block 0 <= b < Q.size()
// A_Y[b], A[b] denote the b-th blocks of A_Y,
// A, resp.

// Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product

// for each b, L.blocks[b], Q[b], Work[b], and A_Y.blocks[b]
// must have the structure described above for `tensorTransposeCongruence'

// TODO: rename this to compute_A_Y, since this Q is
// different from the big Q that gets inverted.

void compute_A_Y(
  const Block_Info &block_info, const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_Y)
{
  A_Y[0].resize(bases_blocks.size() / 2);
  A_Y[1].resize(bases_blocks.size() / 2);

  for(size_t index(0); index < bases_blocks.size(); ++index)
    {
      auto &block(bases_blocks[index]);
      auto &Y_block(Y.blocks[index]);

      El::DistMatrix<El::BigFloat> Y_Q(block),
        A_Y_matrix_temp(block.Width(), block.Width(), block.Grid());
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
           Y_block, block, El::BigFloat(0), Y_Q);

      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), block, Y_Q, El::BigFloat(0), A_Y_matrix_temp);

      const size_t block_size(
        block_info.num_points.at(block_info.block_indices.at(index / 2))),
        dim(block_info.dimensions.at(block_info.block_indices.at(index / 2))),
        parity(index % 2), Q_index(index / 2);
      auto &A_Y_block(A_Y[parity][Q_index]);
      A_Y_block.resize(dim);

      El::MakeSymmetric(El::UpperOrLower::LOWER, A_Y_matrix_temp);

      for(size_t column_block = 0; column_block < dim; ++column_block)
        {
          A_Y_block[column_block].clear();
          A_Y_block[column_block].reserve(dim);
          const size_t column_offset(column_block * block_size);
          for(size_t row_block = 0; row_block < dim; ++row_block)
            {
              const size_t row_offset(row_block * block_size);
              El::DistMatrix<El::BigFloat> submatrix(
                El::LockedView(A_Y_matrix_temp, column_offset, row_offset,
                               block_size, block_size));

              A_Y_block[column_block].emplace_back(block_size, block_size,
                                                   block.Grid());

              El::Transpose(submatrix, A_Y_block[column_block].back());
            }
        }
    }
}
