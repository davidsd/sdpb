#include "sdp_solve/Block_Info.hxx"

void set_bilinear_bases_block(
  const El::Matrix<El::BigFloat> &bilinear_base_local,
  El::DistMatrix<El::BigFloat> &bases_block)
{
  El::BigFloat zero(0);
  for(int64_t row = 0; row < bases_block.LocalHeight(); ++row)
    {
      size_t global_row(bases_block.GlobalRow(row)),
        row_block(global_row / bilinear_base_local.Height());

      for(int64_t column = 0; column < bases_block.LocalWidth(); ++column)
        {
          size_t global_column(bases_block.GlobalCol(column)),
            column_block(global_column / bilinear_base_local.Width());
          bases_block.SetLocal(
            row, column,
            row_block != column_block
              ? zero
              : bilinear_base_local(global_row % bilinear_base_local.Height(),
                                    global_column
                                      % bilinear_base_local.Width()));
        }
    }
}

void set_bases_blocks(
  const Block_Info &block_info,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  const El::Grid &grid)
{
  // Precompute this for compute_bilinear_pairings
  auto pairing_sizes(block_info.bilinear_pairing_block_sizes());
  auto psd_sizes(block_info.psd_matrix_block_sizes());

  auto bilinear(bilinear_bases_local.begin());
  for(auto &block_index : block_info.block_indices)
    {
      for(size_t parity(0); parity != 2; ++parity)
        {
          bases_blocks.emplace_back(psd_sizes[2 * block_index + parity],
                                    pairing_sizes[2 * block_index + parity],
                                    grid);
          auto &block(bases_blocks.back());
          set_bilinear_bases_block(*bilinear, block);
          ++bilinear;
        }
    }
}
