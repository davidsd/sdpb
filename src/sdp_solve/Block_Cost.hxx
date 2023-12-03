#pragma once

#include <vector>
#include <cstdio>

struct Block_Cost
{
  size_t cost, index;
  // Make sure that cost>0.  Otherwise, compute_block_grid_mapping can fail.
  Block_Cost(const size_t &Cost, const size_t &Index)
      : cost(Cost == 0 ? 1 : Cost), index(Index)
  {}
  Block_Cost() = delete;
  bool operator<(const Block_Cost &b) const { return cost < b.cost; }
};
