#include "../Block_Matrix.hxx"

void compute_B_pseudoinverse(const std::vector<size_t> &block_offsets,
                             const std::vector<size_t> &block_indices,
                             const El::Grid &grid, const Block_Matrix &B,
                             Block_Matrix &B_pseudoinverse)
{
  // Copy local elements into a local B
  El::Matrix<El::BigFloat> B_local(block_offsets.back(),
                                   B.blocks.at(0).Width());
  El::Zero(B_local);
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      const size_t block_index(block_indices.at(block));
      auto &B_block(B.blocks.at(block));
      for(int64_t local_row(0); local_row != B_block.LocalHeight();
          ++local_row)
        {
          const int64_t global_row(block_offsets.at(block_index)
                                   + B_block.GlobalRow(local_row));
          for(int64_t local_column(0); local_column != B_block.LocalWidth();
              ++local_column)
            {
              const int64_t global_column(B_block.GlobalRow(local_column));
              B_local(global_row, global_column)
                = B_block.GetLocal(local_row, local_column);
            }
        }
    }

  // Add all elements together from all cores.  For a given element,
  // all of the cores should be zero except one.
  El::mpi::AllReduce(B_local.Buffer(), B_local.Height() * B_local.Width(),
                     El::mpi::SUM, El::mpi::COMM_WORLD);

  // Fill distributed B
  El::DistMatrix<El::BigFloat> B_dist_pinv(B_local.Height(), B_local.Width(),
                                           grid);
  for(int64_t local_row(0); local_row != B_dist_pinv.LocalHeight();
      ++local_row)
    {
      const int64_t global_row(B_dist_pinv.GlobalRow(local_row));
      for(int64_t local_column(0); local_column != B_dist_pinv.LocalWidth();
          ++local_column)
        {
          const int64_t global_column(B_dist_pinv.GlobalRow(local_column));
          B_dist_pinv.SetLocal(local_row, local_column,
                               B_local(global_row, global_column));
        }
    }

  // Compute pseudoinverse
  El::Pseudoinverse(B_dist_pinv);
  
  // Copy to all cores
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> B_star_pinv(B_dist_pinv);

  // Copy to B_pseudoinverse
  B_pseudoinverse.blocks.clear();
  B_pseudoinverse.blocks.reserve(block_indices.size());
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      const size_t block_index(block_indices.at(block));
      B_pseudoinverse.blocks.emplace_back(B.blocks.at(block_index).Height(),
                                          B.blocks.at(block_index).Width(),
                                          grid);
      auto &pseudo_inverse_block(B_pseudoinverse.blocks.back());
      for(int64_t local_row(0);
          local_row != pseudo_inverse_block.LocalHeight(); ++local_row)
        {
          const int64_t global_row(
            block_offsets.at(block_index)
            + pseudo_inverse_block.GlobalRow(local_row));
          for(int64_t local_column(0);
              local_column != pseudo_inverse_block.LocalWidth();
              ++local_column)
            {
              const int64_t global_column(
                pseudo_inverse_block.GlobalCol(local_column));
              pseudo_inverse_block.SetLocal(
                pseudo_inverse_block.LocalRow(local_row),
                pseudo_inverse_block.LocalCol(local_column),
                B_star_pinv.GetLocal(global_row, global_column));
            }
        }
    }
}
