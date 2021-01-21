#include "../../../Block_Diagonal_Matrix.hxx"
#include "../../../Block_Info.hxx"

// Q_X_inv_Q = bilinear_base^T X^{-1} bilinear_base for each block

void compute_Q_X_inv_Q(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &Q_X_inv_Q)
{
  Q_X_inv_Q[0].resize(bases_blocks.size());
  Q_X_inv_Q[1].resize(bases_blocks.size());

  for(size_t index(0); index < bases_blocks.size(); ++index)
    {
      auto &block(bases_blocks[index]);
      auto &X_cholesky_block(X_cholesky.blocks[index]);
      El::DistMatrix<El::BigFloat> temp_space(block),
        Q_X_inv_Q_matrix(block.Width(), block.Width(), block.Grid());
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), X_cholesky_block, temp_space);

      // We have to set this to zero because the values can be NaN.
      // Multiplying 0*NaN = NaN.
      Zero(Q_X_inv_Q_matrix);
      Syrk(El::UpperOrLowerNS::LOWER, El::Orientation::TRANSPOSE,
           El::BigFloat(1), temp_space, El::BigFloat(0), Q_X_inv_Q_matrix);
      El::MakeSymmetric(El::UpperOrLower::LOWER, Q_X_inv_Q_matrix);

      const size_t block_size(
        block_info.num_points.at(block_info.block_indices.at(index / 2))),
        dim(block_info.dimensions.at(block_info.block_indices.at(index / 2)));

      const size_t parity(index % 2), Q_index(index / 2);
      auto &Q_X_inv_Q_block(Q_X_inv_Q[parity][Q_index]);
      Q_X_inv_Q_block.resize(dim);
      for(size_t column_block = 0; column_block < dim; ++column_block)
        {
          Q_X_inv_Q_block[column_block].clear();
          Q_X_inv_Q_block[column_block].reserve(dim);
          const size_t column_offset(column_block * block_size);
          for(size_t row_block = 0; row_block < dim; ++row_block)
            {
              const size_t row_offset(row_block * block_size);
              El::DistMatrix<El::BigFloat> submatrix(
                El::View(Q_X_inv_Q_matrix, column_offset, row_offset,
                         block_size, block_size));

              Q_X_inv_Q_block[column_block].emplace_back(
                block_size, block_size, block.Grid());
              El::Copy(submatrix, Q_X_inv_Q_block[column_block].back());
            }
        }
    }
}
