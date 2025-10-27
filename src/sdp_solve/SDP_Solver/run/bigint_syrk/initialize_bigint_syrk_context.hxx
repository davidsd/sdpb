#pragma once

#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/SDP.hxx"

inline BigInt_Shared_Memory_Syrk_Context
initialize_bigint_syrk_context(const Block_Info &block_info, const SDP &sdp,
                               const size_t max_shared_memory_bytes,
                               const Verbosity verbosity)
{
  const auto group_index = block_info.node_group_index();
  const auto precision = El::gmp::Precision();

  const auto block_width = sdp.dual_objective_b.Height(); // = N

  std::vector<int> blocks_height_per_group;
  for(const auto &group :
      block_info.block_mapping.mapping.at(block_info.node_index))
    {
      size_t height = 0;
      for(const auto &block_index : group.block_indices)
        {
          // P' = height of the current block of sdp.free_var_matrix (aka B)
          height += block_info.get_schur_block_size(block_index);
        }
      blocks_height_per_group.push_back(height);
    }

  return BigInt_Shared_Memory_Syrk_Context(
    block_info.node_comm, group_index, precision, max_shared_memory_bytes,
    blocks_height_per_group, block_width, block_info.block_indices, verbosity);
}
