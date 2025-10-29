#include "Block_Mapping.hxx"

Block_Mapping::Block_Mapping(const std::vector<std::vector<Block_Map>> &mapping)
    : mapping(mapping)
{
  const auto num_nodes = mapping.size();
  node_rank_to_group.resize(num_nodes);
  block_index_node_to_global.resize(num_nodes);

  for(size_t node_index = 0; node_index < mapping.size(); ++node_index)
    {
      auto &rank_to_group = node_rank_to_group.at(node_index);
      auto &block_indices = block_index_node_to_global.at(node_index);
      for(size_t group_index = 0; group_index < mapping.at(node_index).size();
          ++group_index)
        {
          const auto &group = mapping.at(node_index).at(group_index);
          rank_to_group.insert(rank_to_group.end(), group.num_procs,
                               group_index);
          block_indices.insert(block_indices.end(),
                               group.block_indices.begin(),
                               group.block_indices.end());
        }
    }
}
