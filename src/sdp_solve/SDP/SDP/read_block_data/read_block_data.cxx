#include "SDP_Block_Data.hxx"
#include "pmp2sdp/write_sdp.hxx"
#include "sdp_solve/Archive_Reader.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/copy_matrix.hxx"

#include <unordered_map>

namespace fs = std::filesystem;

namespace
{
  Block_File_Format get_block_format(const fs::path &block_path)
  {
    const auto &extension = block_path.extension();
    if(extension == ".json")
      return Block_File_Format::json;
    if(extension == ".bin")
      return Block_File_Format::bin;
    RUNTIME_ERROR("Unknown block file extension: ", block_path);
  }

  // sdp_block_local in initialized only at comm.Rank() == 0
  // Data from sdp_block_local is sent to DistMatrices in sdp
  void set_sdp_from_root(const El::Grid &grid, const Block_Info &block_info,
                         const SDP_Block_Data &sdp_block_local, SDP &sdp)
  {
    const auto &comm = grid.Comm();

    int index = sdp_block_local.block_index_local;
    El::mpi::Broadcast(index, 0, comm);
    if(index == -1)
      RUNTIME_ERROR("block_data file not found for one of the indices");

    const size_t block_index = block_info.block_indices.at(index);

    // sdp.primal_objective_c
    {
      auto &c(sdp.primal_objective_c.blocks.at(index));
      c.SetGrid(grid);
      c.Resize(block_info.get_schur_block_size(block_index), 1);
      copy_matrix_from_root(sdp_block_local.primal_objective_c, c, comm);
    }

    // sdp.free_var_matrix
    {
      auto &B(sdp.free_var_matrix.blocks.at(index));
      B.SetGrid(grid);
      // B block has size P'*N
      // NB: sdp.dual_objective_b must be initialized at this moment!
      // This is done in practice in SDP constructor, but no guaranteed generally.
      // TODO initialize it in Block_Info, for consistency?
      B.Resize(block_info.get_schur_block_size(block_index),
               sdp.dual_objective_b.Height());
      copy_matrix_from_root(sdp_block_local.constraint_matrix, B, comm);
    }

    // sdp.bilinear_bases and sdp.bases_blocks
    for(const size_t parity : {0, 1})
      {
        const size_t bilinear_index_local = 2 * index + parity;

        // Set sdp.bilinear_bases:
        {
          const auto &bilinear_bases_local
            = sdp_block_local.bilinear_bases[parity];
          auto &bilinear_bases = sdp.bilinear_bases.at(bilinear_index_local);

          bilinear_bases.SetGrid(grid);
          const auto height
            = block_info.get_bilinear_bases_height(block_index, parity);
          const auto width
            = block_info.get_bilinear_bases_width(block_index, parity);
          bilinear_bases.Resize(height, width);
          copy_matrix_from_root(bilinear_bases_local, bilinear_bases, comm);
        }

        // Set sdp.bases_blocks:
        {
          auto &bases_block_local = sdp_block_local.bases_blocks[parity];

          const auto height
            = block_info.get_psd_matrix_block_size(block_index, parity);
          const auto width
            = block_info.get_bilinear_pairing_block_size(block_index, parity);

          auto &bases_block = sdp.bases_blocks.at(bilinear_index_local);
          bases_block.SetGrid(grid);
          bases_block.Resize(height, width);
          copy_matrix_from_root(bases_block_local, bases_block, comm);
        }
      }
  }
}

void read_block_data(const fs::path &sdp_path, const El::Grid &grid,
                     const Block_Info &block_info, SDP &sdp, Timers &timers)
{
  Scoped_Timer timer(timers, "read_block_data");

  ASSERT(fs::exists(sdp_path), "SDP path does not exist: ", sdp_path);

  const size_t num_blocks(block_info.block_indices.size());
  sdp.primal_objective_c.blocks.resize(num_blocks);
  sdp.free_var_matrix.blocks.resize(num_blocks);
  sdp.bilinear_bases.resize(2 * num_blocks);
  sdp.bases_blocks.resize(2 * num_blocks);

  const auto &comm = grid.Comm();

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

      std::unique_ptr<Archive_Reader> reader;
      if(comm.Rank() == 0)
        reader = std::make_unique<Archive_Reader>(sdp_path);

      for(size_t processed_count = 0; processed_count < num_blocks;
          ++processed_count)
        {
          Scoped_Timer count_timer(timers, std::to_string(processed_count));

          // Initialized on root only
          SDP_Block_Data sdp_block_local;

          if(comm.Rank() == 0)
            {
              // Find and read next entry with one of the block indices required
              while(reader->next_entry())
                {
                  const fs::path curr_block_path
                    = archive_entry_pathname(reader->entry_ptr);
                  int index = get_index(curr_block_path);
                  if(index == -1)
                    continue;

                  Scoped_Timer parse_timer(timers, "parse");
                  Block_File_Format format = get_block_format(curr_block_path);
                  std::istream stream(reader.get());
                  sdp_block_local
                    = SDP_Block_Data(stream, format, index, block_info);
                  break;
                }
            }

          Scoped_Timer sync_timer(timers, "synchronize");
          set_sdp_from_root(grid, block_info, sdp_block_local, sdp);
        }
    }
  else
    {
      for(size_t index(0); index != num_blocks; ++index)
        {
          Scoped_Timer count_timer(timers, std::to_string(index));

          // Initialized on root only
          SDP_Block_Data sdp_block_local;

          if(comm.Rank() == 0)
            {
              Scoped_Timer parse_timer(timers, "parse");
              fs::path block_path(
                sdp_path
                / ("block_data_"
                   + std::to_string(block_info.block_indices.at(index))
                   + ".bin"));
              if(!exists(block_path))
                block_path.replace_extension(".json");
              ASSERT(exists(block_path), "Block not found: ", block_path);

              Block_File_Format format = get_block_format(block_path);
              std::ifstream block_stream(block_path, std::ios::binary);
              sdp_block_local
                = SDP_Block_Data(block_stream, format, index, block_info);
            }

          Scoped_Timer sync_timer(timers, "synchronize");
          set_sdp_from_root(grid, block_info, sdp_block_local, sdp);
        }
    }
}
