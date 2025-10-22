#pragma once

#include "Block_Map.hxx"

#include <vector>

struct Block_Mapping
{
  // Node -> Group -> Block_Map(num_procs, block_indices)
  std::vector<std::vector<Block_Map>> mapping;
  // Node -> Rank (on a node) -> Group
  std::vector<std::vector<size_t>> node_rank_to_group;

  explicit Block_Mapping(const std::vector<std::vector<Block_Map>> &mapping)
      : mapping(mapping), node_rank_to_group(to_group_indices(mapping))
  {}

private:
  static std::vector<std::vector<size_t>>
  to_group_indices(const std::vector<std::vector<Block_Map>> &mapping)
  {
    std::vector<std::vector<size_t>> result(mapping.size());
    for(size_t node_index = 0; node_index < mapping.size(); ++node_index)
      {
        for(size_t group_index = 0;
            group_index < mapping.at(node_index).size(); ++group_index)
          {
            const auto &group = mapping.at(node_index).at(group_index);
            result.at(node_index)
              .insert(result.at(node_index).end(), group.num_procs,
                      group_index);
          }
      }
    return result;
  }
};
