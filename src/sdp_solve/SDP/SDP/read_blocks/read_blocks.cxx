#include "../set_bases_blocks.hxx"
#include "../../../SDP.hxx"
#include "../../../Archive_Reader.hxx"

#include <algorithm>

namespace
{
  struct Entry
  {
    std::string s;
    size_t index;
    Entry(const std::string &S, const size_t &Index): s(S), index(Index) {}
    bool operator<(const Entry &b) const
    {
      return s < b.s;
    }
  };
  inline bool operator<(const std::string &a, const Entry &b) { return a < b.s; }
  inline bool operator<(const Entry &a, const std::string &b) { return a.s < b; }
}

void read_block_stream(
  const El::Grid &grid, const size_t &index, std::istream &block_stream,
  SDP &sdp, std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local);

void read_blocks(const boost::filesystem::path &sdp_path, const El::Grid &grid,
                 const Block_Info &block_info, SDP &sdp)
{
  const size_t num_blocks(block_info.block_indices.size());
  sdp.primal_objective_c.blocks.resize(num_blocks);
  sdp.free_var_matrix.blocks.resize(num_blocks);

  std::vector<El::Matrix<El::BigFloat>> bilinear_bases_local(2 * num_blocks);

  if(boost::filesystem::is_regular_file(sdp_path))
    {
      // The archive could have 10's of thousands of entries, so we
      // make a fast lookup table.
      std::vector<Entry> block_names;
      block_names.reserve(block_info.block_indices.size());
      for(size_t index(0); index != num_blocks; ++index)
        {
          block_names.emplace_back(
            "block_" + std::to_string(block_info.block_indices[index])
              + ".json",
            index);
        }
      std::sort(block_names.begin(), block_names.end());
      
      Archive_Reader reader(sdp_path);
      while(reader.next_entry())
        {
          const std::string entry(archive_entry_pathname(reader.entry_ptr));
          auto block_entry(std::equal_range(block_names.begin(),
                                            block_names.end(), entry));
          if(block_entry.first != block_entry.second)
            {
              std::istream stream(&reader);
              read_block_stream(grid, block_entry.first->index, stream, sdp,
                                bilinear_bases_local);
            }
        }
    }
  else
    {
      for(size_t index(0); index != num_blocks; ++index)
        {
          const boost::filesystem::path block_path(
            sdp_path
            / ("block_" + std::to_string(block_info.block_indices.at(index))
               + ".json"));
          boost::filesystem::ifstream block_stream(block_path);
          read_block_stream(grid, index, block_stream, sdp,
                            bilinear_bases_local);
        }
    }

  sdp.bilinear_bases.reserve(num_blocks);
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
