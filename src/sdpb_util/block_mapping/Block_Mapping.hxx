#pragma once

#include "Block_Map.hxx"

#include <map>
#include <ostream>
#include <vector>

struct Block_Location
{
  size_t node_index;
  // MPI group on a node
  size_t group_index;
  size_t block_index_local;
  size_t block_index_node;
  size_t block_index_global;

  Block_Location() = delete;
  Block_Location(size_t node_index, size_t group_index,
                 size_t block_index_local, size_t block_index_node,
                 size_t block_index_global);
  friend bool operator==(const Block_Location &lhs, const Block_Location &rhs);
  friend bool operator!=(const Block_Location &lhs, const Block_Location &rhs);
  friend std::ostream &operator<<(std::ostream &os, const Block_Location &loc);
};

struct Block_Mapping
{
  // Node -> Group -> Block_Map(num_procs, block_indices)
  std::vector<std::vector<Block_Map>> mapping;
  // block index (global) -> location
  std::vector<Block_Location> block_locations;
  // Node -> block locations
  std::vector<std::vector<Block_Location>> node_block_locations;
  // Node -> Group -> block locations
  std::vector<std::vector<std::vector<Block_Location>>>
    node_group_block_locations;
  // Node -> Rank (on a node) -> Group
  std::vector<std::vector<size_t>> node_rank_to_group;
  explicit Block_Mapping(const std::vector<std::vector<Block_Map>> &mapping);

  size_t num_nodes() const;
  size_t num_groups(size_t node_index) const;

private:
  void validate() const;
};
