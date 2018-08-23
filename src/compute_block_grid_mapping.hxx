#pragma once

#include "Block_Cost.hxx"
#include "Block_Map.hxx"

std::vector<std::vector<Block_Map>>
compute_block_grid_mapping(const size_t &procs_per_node,
                           const size_t &num_nodes,
                           const std::vector<Block_Cost> &block_costs);
