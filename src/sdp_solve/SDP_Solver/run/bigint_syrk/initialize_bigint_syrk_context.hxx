#pragma once

#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/Block_Matrix.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdpb_util/assert.hxx"

// parent_comm is split into child_comms.
// We assign child_index = [0..num_children) to each child_comm
// and return {child_index, [child_comm.Size() for each child], num_children}
inline std::tuple<size_t, std::vector<int>, size_t>
get_child_index_sizes_and_count(const El::mpi::Comm &parent_comm,
                                const El::mpi::Comm &child_comm)
{
  ASSERT(parent_comm.Size() > 0);
  ASSERT(child_comm.Size() > 0);

  size_t child_index = std::numeric_limits<size_t>::max();
  size_t num_children = 0;
  std::vector<int> child_comm_sizes;

  size_t curr_index = 0;
  for(int curr_rank = 0; curr_rank < parent_comm.Size(); ++curr_rank)
    {
      // When we are on the first rank of the parent, set index
      if(curr_rank == parent_comm.Rank() && child_comm.Rank() == 0)
        {
          child_index = curr_index;
          curr_index++;
          child_comm_sizes.push_back(child_comm.Size());
        }
      // All ranks of a child_comm have the same index
      El::mpi::Broadcast(child_index, 0, child_comm);
      // Synchronize curr_index
      El::mpi::Broadcast(curr_index, curr_rank, parent_comm);
      // Synchronize child_comm_sizes
      auto curr_child_index = curr_index - 1;
      child_comm_sizes.resize(curr_child_index + 1, 0);
      El::mpi::Broadcast(child_comm_sizes.at(curr_child_index), curr_rank,
                         parent_comm);
    }
  num_children = curr_index;

  ASSERT_EQUAL(child_comm_sizes.size(), num_children);
  ASSERT(num_children > 0);
  ASSERT(child_index < num_children, DEBUG_STRING(child_index),
         DEBUG_STRING(num_children));
  ASSERT_EQUAL(parent_comm.Size(),
               std::accumulate(child_comm_sizes.begin(),
                               child_comm_sizes.end(), (size_t)0));

  return {child_index, child_comm_sizes, num_children};
}

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
  Grouped_Block_Size_Info(const Environment &env, const Block_Info &block_info,
                          const SDP &sdp)
      : node_comm(env.comm_shared_mem), group_comm(block_info.mpi_comm.value)
  {
    block_width = sdp.dual_objective_b.Height(); // = N

    std::tie(group_index, group_comm_sizes, num_groups)
      = get_child_index_sizes_and_count(node_comm, group_comm);
    blocks_height_per_group.resize(num_groups, 0);

    for(size_t i = 0; i < num_groups; ++i)
      {
        if(i == group_index && group_comm.Rank() == 0)
          {
            // Total height for all blocks on group_comm.
            // Since all ranks of group_comm share the same blocks,
            // It's enough to calculate it on group root only.
            blocks_height_per_group.at(i) = std::transform_reduce(
              sdp.free_var_matrix.blocks.begin(),
              sdp.free_var_matrix.blocks.end(), 0, std::plus{},
              [](const auto &block) { return block.Height(); });
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
