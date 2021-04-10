#include "../set_bases_blocks.hxx"
#include "../../../SDP.hxx"
#include "../../../Archive_Reader.hxx"

void read_block_stream(
  const El::Grid &grid, std::istream &block_stream, SDP &sdp,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local);

void read_blocks(const boost::filesystem::path &sdp_path, const El::Grid &grid,
                 const Block_Info &block_info, SDP &sdp)
{
  sdp.primal_objective_c.blocks.reserve(block_info.block_indices.size());

  std::vector<El::Matrix<El::BigFloat>> bilinear_bases_local;
  bilinear_bases_local.reserve(2 * block_info.block_indices.size());
  sdp.free_var_matrix.blocks.reserve(block_info.block_indices.size());

  for(auto &block_index : block_info.block_indices)
    {
      if(boost::filesystem::is_regular_file(sdp_path))
        {
          // TODO: This is going to reopen the zip file many, many
          // times.
          const std::string block_name("block_" + std::to_string(block_index)
                                       + ".json");
          Archive_Reader reader(sdp_path);
          while(reader.next_entry())
            {
              if(block_name == archive_entry_pathname(reader.entry_ptr))
                {
                  std::istream stream(&reader);
                  read_block_stream(grid, stream, sdp, bilinear_bases_local);
                }
            }
        }
      else
        {
          const boost::filesystem::path block_path(
            sdp_path / ("block_" + std::to_string(block_index) + ".json"));
          boost::filesystem::ifstream block_stream(block_path);
          read_block_stream(grid, block_stream, sdp, bilinear_bases_local);
        }
    }

  sdp.bilinear_bases.reserve(bilinear_bases_local.size());
  for(auto &local : bilinear_bases_local)
    {
      sdp.bilinear_bases.emplace_back(local.Height(), local.Width(), grid);
      auto &dist(sdp.bilinear_bases.back());
      for(int64_t local_row(0); local_row < dist.LocalHeight(); ++local_row)
        {
          El::Int global_row(dist.GlobalRow(local_row));
          for(int64_t local_column(0); local_column < dist.LocalWidth();
              ++local_column)
            {
              El::Int global_column(dist.GlobalCol(local_column));
              dist.SetLocal(local_row, local_column,
                            local(global_row, global_column));
            }
        }
    }
  set_bases_blocks(block_info, bilinear_bases_local, sdp.bases_blocks, grid);
}
