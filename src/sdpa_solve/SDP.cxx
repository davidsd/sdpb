#include "SDP.hxx"

#include "parse_sdpa/parse_sdpa.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/copy_matrix.hxx"

#include <filesystem>

namespace fs = std::filesystem;

namespace Sdpb::Sdpa
{
  SDP::SDP(const fs::path &sdp_path, const Block_Info &block_info,
           const El::Grid &grid, Timers &timers)
      : sdp_block_F_0(block_info.block_dimensions, block_info.block_indices,
                      grid)
  {
    Scoped_Timer timer(timers, "SDP_ctor");

    const auto m = block_info.primal_dimension;

    const auto &comm = block_info.mpi_comm.value;

    // For each MPI group, rank=0 parses objective_c and the blocks for F_0..F_m that belong to this group.
    // Other ranks do not parse anything.
    const auto parse_result = [&] {
      if(comm.Rank() == 0)
        {
          Scoped_Timer parse_timer(timers, "parse");
          std::set block_indices_to_parse(block_info.block_indices.begin(),
                                          block_info.block_indices.end());
          const auto should_parse_block = [&](const size_t index) {
            return block_indices_to_parse.find(index)
                   != block_indices_to_parse.end();
          };
          return read_sdpa(sdp_path, should_parse_block);
        }
      return SDPA_File_Parse_Result();
    }();

    {
      Scoped_Timer sync_timer(timers, "synchronize");

      // primal_objective_c
      {
        // All elements are stored at rank=0
        auto c_circ_circ
          = El::DistMatrix<El::BigFloat, El::CIRC, El::CIRC>(m, 1);
        if(El::mpi::Rank() == 0)
          {
            for(size_t i = 0; i < m; ++i)
              c_circ_circ.Set(i, 0, parse_result.c_objective.at(i));
          }
        // Copy all elements from rank=0 to all ranks.
        primal_objective_c = c_circ_circ;
      }

      // Prepare F_0 and F_i structure
      {
        sdp_block_F_0 = Block_Diagonal_Matrix(block_info.block_dimensions,
                                              block_info.block_indices, grid);
        sdp_blocks_F.reserve(m);
        for(size_t i = 0; i < m; ++i)
          {
            sdp_blocks_F.emplace_back(block_info.block_dimensions,
                                      block_info.block_indices, grid);
          }
      }

      // Synchronize F_0 and F_i blocks
      for(size_t local_block_index = 0;
          local_block_index < block_info.num_blocks_local();
          ++local_block_index)
        {
          const auto block_index
            = block_info.block_indices.at(local_block_index);
          // Run over F_0, F_1..F_m
          // NB: sdp_blocks_F[i] corresponds to parsed_block[i+1]
          // since F_0 is not in sdp_blocks_F
          for(size_t i = 0; i <= m; ++i)
            {
              const auto &source
                = comm.Rank() == 0
                    ? parse_result.parsed_blocks.at(block_index).at(i)
                    : El::Matrix<El::BigFloat>();
              auto &F = i == 0 ? sdp_block_F_0 : sdp_blocks_F.at(i - 1);
              auto &dest = F.blocks.at(local_block_index);
              copy_matrix_from_root(source, dest, comm);
            }
        }
    }

    Scoped_Timer validate_timer(timers, "validate");
    validate(block_info);
  }

  size_t SDP::primal_dimension() const
  {
    return sdp_blocks_F.size();
  }

  void SDP::validate(const Block_Info &block_info) const noexcept(false)
  {
    const auto error_prefix = "Invalid SDP: ";

    // Check primal_objective_c
    {
      ASSERT_EQUAL(primal_objective_c.Height(), this->sdp_blocks_F.size(),
                   error_prefix,
                   "c_1..c_m and F_1..Fm should have the same length");
      ASSERT_EQUAL(primal_objective_c.Width(), 1, error_prefix,
                   "c should be a vector");
    }

    for(size_t local_block_index = 0;
        local_block_index < block_info.num_blocks_local(); ++local_block_index)
      {
        const auto block_index
          = block_info.block_indices.at(local_block_index);

        const auto error_prefix_index
          = El::BuildString(error_prefix, "block_index = ", block_index, ": ");

        const auto &comm = block_info.mpi_comm.value;
        const auto dim = block_info.block_dimensions.at(block_index);

        {
          const auto &block = sdp_block_F_0.blocks.at(local_block_index);
          ASSERT(El::mpi::Congruent(block.DistComm(), comm),
                 error_prefix_index, "wrong MPI communicator for F_0");
          ASSERT_EQUAL(block.Height(), dim, error_prefix_index);
          ASSERT_EQUAL(block.Width(), dim, error_prefix_index);
        }

        for(const auto &F : sdp_blocks_F)
          {
            const auto &block = F.blocks.at(local_block_index);
            ASSERT(El::mpi::Congruent(block.DistComm(), comm),
                   error_prefix_index,
                   "wrong MPI communicator for sdp_blocks_F");
            ASSERT_EQUAL(block.Height(), dim, error_prefix_index);
            ASSERT_EQUAL(block.Width(), dim, error_prefix_index);
          }
      }
  }
}
