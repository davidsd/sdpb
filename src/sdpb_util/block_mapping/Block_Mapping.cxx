#include "Block_Mapping.hxx"

#include "sdpb_util/assert.hxx"

Block_Location::Block_Location(const size_t node_index,
                               const size_t group_index,
                               const size_t block_index_local,
                               const size_t block_index_node,
                               const size_t block_index_global)
    : node_index(node_index),
      group_index(group_index),
      block_index_local(block_index_local),
      block_index_node(block_index_node),
      block_index_global(block_index_global)
{}
bool operator==(const Block_Location &lhs, const Block_Location &rhs)
{
  return lhs.node_index == rhs.node_index && lhs.group_index == rhs.group_index
         && lhs.block_index_local == rhs.block_index_local
         && lhs.block_index_node == rhs.block_index_node
         && lhs.block_index_global == rhs.block_index_global;
}
bool operator!=(const Block_Location &lhs, const Block_Location &rhs)
{
  return !(lhs == rhs);
}
std::ostream &operator<<(std::ostream &os, const Block_Location &loc)
{
  return os << "node_index: " << loc.node_index
            << " group_index: " << loc.group_index
            << " block_index_local: " << loc.block_index_local
            << " block_index_node: " << loc.block_index_node
            << " block_index_global: " << loc.block_index_global;
}

Block_Mapping::Block_Mapping(const std::vector<std::vector<Block_Map>> &mapping)
    : mapping(mapping)
{
  const auto num_nodes = this->num_nodes();
  node_rank_to_group.resize(num_nodes);

  // Global block index -> location
  std::map<size_t, Block_Location> block_location_map;

  node_block_locations.resize(num_nodes);
  node_group_block_locations.resize(num_nodes);

  size_t num_blocks = 0;
  for(size_t node_index = 0; node_index < num_nodes; ++node_index)
    {
      auto &rank_to_group = node_rank_to_group.at(node_index);
      const auto num_groups = this->num_groups(node_index);
      node_group_block_locations.at(node_index).resize(num_groups);

      size_t block_index_node = 0;
      for(size_t group_index = 0; group_index < num_groups; ++group_index)
        {
          const auto &group = mapping.at(node_index).at(group_index);
          rank_to_group.insert(rank_to_group.end(), group.num_procs,
                               group_index);
          for(size_t block_index_local = 0;
              block_index_local < group.block_indices.size();
              ++block_index_local, ++block_index_node, ++num_blocks)
            {
              const size_t block_index_global
                = group.block_indices.at(block_index_local);
              Block_Location block_location(
                node_index, group_index, block_index_local, block_index_node,
                block_index_global);
              auto [it, inserted] = block_location_map.insert(
                {block_index_global, block_location});
              ASSERT(inserted, "Duplicate block_index=", block_index_global,
                     DEBUG_STRING(node_index), DEBUG_STRING(group_index));
              node_block_locations.at(node_index).push_back(block_location);
              node_group_block_locations.at(node_index)
                .at(group_index)
                .push_back(block_location);
            }
        }
    }

  block_locations.reserve(num_blocks);
  for(size_t global_block_index = 0; global_block_index < num_blocks;
      ++global_block_index)
    {
      auto it = block_location_map.find(global_block_index);
      ASSERT(it != block_location_map.end(),
             "Block not found! index=", global_block_index);
      block_locations.push_back(it->second);
    }

  validate();
}
size_t Block_Mapping::num_nodes() const
{
  return mapping.size();
}
size_t Block_Mapping::num_groups(const size_t node_index) const
{
  return mapping.at(node_index).size();
}
void Block_Mapping::validate() const
{
  ASSERT_EQUAL(mapping.size(), num_nodes());
  ASSERT_EQUAL(node_block_locations.size(), num_nodes());
  ASSERT_EQUAL(node_group_block_locations.size(), num_nodes());
  for(size_t node_index = 0; node_index < num_nodes(); ++node_index)
    {
      size_t block_index_node = 0;
      ASSERT_EQUAL(mapping.at(node_index).size(), num_groups(node_index));
      ASSERT_EQUAL(node_group_block_locations.at(node_index).size(),
                   num_groups(node_index));
      for(size_t group_index = 0; group_index < num_groups(node_index);
          ++group_index)
        {
          const auto &block_map = mapping.at(node_index).at(group_index);
          ASSERT_EQUAL(
            block_map.block_indices.size(),
            node_group_block_locations.at(node_index).at(group_index).size());
          for(size_t block_index_local = 0;
              block_index_local < block_map.block_indices.size();
              ++block_index_local, ++block_index_node)
            {
              const auto block_index_global
                = block_map.block_indices.at(block_index_local);
              Block_Location location(node_index, group_index,
                                      block_index_local, block_index_node,
                                      block_index_global);

              ASSERT_EQUAL(location, block_locations.at(block_index_global));
              ASSERT_EQUAL(
                location,
                node_block_locations.at(node_index).at(block_index_node));
              ASSERT_EQUAL(location, node_group_block_locations.at(node_index)
                                       .at(group_index)
                                       .at(block_index_local));
            }
        }
    }
  // TODO check also node_rank_to_group
}
