#pragma once

#include "Block_Map.hxx"

#include <vector>

struct Block_Mapping
{
  // Node -> Group -> Block_Map(num_procs, block_indices)
  std::vector<std::vector<Block_Map>> mapping;
  // Node -> Rank (on a node) -> Group
  std::vector<std::vector<size_t>> node_rank_to_group;
  // node -> block_index(on node) -> block index
  std::vector<std::vector<size_t>> block_index_node_to_global;

  explicit Block_Mapping(const std::vector<std::vector<Block_Map>> &mapping);
};
