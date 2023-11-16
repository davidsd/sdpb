#include "../set_bases_blocks.hxx"
#include "../../../SDP.hxx"
#include "../../../Archive_Reader.hxx"
#include "../../../../sdp_convert.hxx"

#include <algorithm>
#include <unordered_map>

namespace fs = std::filesystem;

namespace
{
  Block_File_Format get_block_format(const fs::path &block_path)
  {
    auto extension = block_path.extension();
    if(extension == ".json")
      return Block_File_Format::json;
    else if(extension == ".bin")
      return Block_File_Format::bin;
    El::RuntimeError("Unknown block file extension: ", block_path);
  }
}

void read_block_stream(
  const El::Grid &grid, const size_t &index, std::istream &block_stream,
  Block_File_Format format, SDP &sdp,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local);

void read_blocks(const fs::path &sdp_path, const El::Grid &grid,
                 const Block_Info &block_info, SDP &sdp)
{
  if(!fs::exists(sdp_path))
    {
      El::RuntimeError("SDP path '" + sdp_path.string() + "' does not exist");
    }

  const size_t num_blocks(block_info.block_indices.size());
  sdp.primal_objective_c.blocks.resize(num_blocks);
  sdp.free_var_matrix.blocks.resize(num_blocks);

  std::vector<El::Matrix<El::BigFloat>> bilinear_bases_local(2 * num_blocks);

  if(fs::is_regular_file(sdp_path))
    {
      // The archive could have 10's of thousands of entries, so we
      // make a fast lookup table.
      std::unordered_map<std::string, size_t> index_by_blockname;
      for(size_t index(0); index != num_blocks; ++index)
        {
          index_by_blockname.emplace(
            "block_data_" + std::to_string(block_info.block_indices[index]),
            index);
        }
      auto get_index
        = [&index_by_blockname](const fs::path &block_path) -> int {
        auto name = fs::path(block_path).replace_extension().string();
        auto index_it = index_by_blockname.find(name);
        if(index_it == index_by_blockname.end())
          return -1;
        return index_it->second;
      };

      Archive_Reader reader(sdp_path);
      size_t processed_count = 0;
      while(reader.next_entry())
        {
          const fs::path curr_block_path
            = archive_entry_pathname(reader.entry_ptr);
          int index = get_index(curr_block_path);
          if(index == -1)
            continue;
          Block_File_Format format = get_block_format(curr_block_path);
          std::istream stream(&reader);
          read_block_stream(grid, index, stream, format, sdp,
                            bilinear_bases_local);
          processed_count++;
          if(processed_count == num_blocks)
            break;
        }
    }
  else
    {
      for(size_t index(0); index != num_blocks; ++index)
        {
          fs::path block_path(
            sdp_path
            / ("block_data_"
               + std::to_string(block_info.block_indices.at(index)) + ".bin"));
          if(!exists(block_path))
            block_path.replace_extension(".json");
          if(!exists(block_path))
            El::RuntimeError("Block not found: ", block_path);

          Block_File_Format format = get_block_format(block_path);
          std::ifstream block_stream(block_path, std::ios::binary);
          read_block_stream(grid, index, block_stream, format, sdp,
                            bilinear_bases_local);
        }
    }

  sdp.bilinear_bases.reserve(num_blocks);
  for(auto &local : bilinear_bases_local)
    {
      sdp.bilinear_bases.emplace_back(local.Height(), local.Width(), grid);
      auto &dist(sdp.bilinear_bases.back());
      for(El::Int local_row(0); local_row < dist.LocalHeight(); ++local_row)
        {
          El::Int global_row(dist.GlobalRow(local_row));
          for(El::Int local_column(0); local_column < dist.LocalWidth();
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
