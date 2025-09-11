#pragma once

#include "sdp_solve/Block_Matrix/Block_Matrix.hxx"
#include "sdpa_solve/Block_Info.hxx"
#include "sdpa_solve/SDP.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/initialize_bigint_syrk_context.hxx"
#include "sdpb_util/assert.hxx"

namespace Sdpb::Sdpa
{
  struct Grouped_Block_Size_Info
  {
    // All ranks on a node
    El::mpi::Comm node_comm;
    // node_comm is split into groups, each SDP block is assigned to one group
    El::mpi::Comm group_comm;
    // current group index in node, [0, num_groups)
    size_t group_index;
    // Number of groups on the node
    size_t num_groups;
    // Size of group_comm for each group
    std::vector<int> group_comm_sizes;

    size_t block_width;
    // Total block height for blocks owned by a given group_comm
    std::vector<El::Int> blocks_height_per_group;

    [[nodiscard]]
    Grouped_Block_Size_Info(const Environment &env,
                            const Block_Info &block_info, const SDP &sdp)
        : node_comm(env.comm_shared_mem), group_comm(block_info.mpi_comm.value)
    {
      // The input matrix is P = (vec(G_1), vec(G_2),...vec(G_m))
      // Where G_p is a block-diagonal matrix of the same size as F_i.
      // This means that width(P) = primal_dimension
      // and block_heights = [(dim_i * dim_i) for i in SDP blocks]
      block_width = sdp.primal_dimension();

      std::tie(group_index, group_comm_sizes, num_groups)
        = get_child_index_sizes_and_count(node_comm, group_comm);
      blocks_height_per_group.resize(num_groups, 0);

      for(size_t i = 0; i < num_groups; ++i)
        {
          if(i == group_index && group_comm.Rank() == 0)
            {
              const auto &F_0 = sdp.sdp_block_F_0;
              // Total height for all blocks on group_comm.
              // Since all ranks of group_comm share the same blocks,
              // It's enough to calculate it on group root only.
              blocks_height_per_group.at(i) = std::transform_reduce(
                F_0.blocks.begin(), F_0.blocks.end(), 0, std::plus{},
                [](const auto &block) {
                  return block.Height() * block.Width();
                });
            }
        }
      El::mpi::AllReduce(blocks_height_per_group.data(), num_groups,
                         El::mpi::SUM, node_comm);
    }
  };

  inline BigInt_Shared_Memory_Syrk_Context
  initialize_bigint_syrk_context(const Environment &env,
                                 const Block_Info &block_info, const SDP &sdp,
                                 const size_t max_shared_memory_bytes,
                                 const Verbosity verbosity)
  {
    const Grouped_Block_Size_Info info(env, block_info, sdp);

    return BigInt_Shared_Memory_Syrk_Context(
      env.comm_shared_mem, info.group_index, info.group_comm_sizes,
      El::gmp::Precision(), max_shared_memory_bytes,
      info.blocks_height_per_group, info.block_width, block_info.block_indices,
      verbosity);
  }
}
