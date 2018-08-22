#pragma once

#include "Block_Cost.hxx"
#include "Block_Map.hxx"

std::vector<Block_Map>
compute_block_grid_mapping(const size_t &num_procs,
                           const std::vector<Block_Cost> &costs);
