#pragma once

#include <cstddef>
#include <tuple>

struct Block_Cost
{
  size_t cost, index;
  // Make sure that cost>0.  Otherwise, compute_block_grid_mapping can fail.
  Block_Cost(const size_t &Cost, const size_t &Index)
      : cost(Cost == 0 ? 1 : Cost), index(Index)
  {}
  Block_Cost() = delete;
  bool operator<(const Block_Cost &b) const
  {
    return std::tie(cost, index) < std::tie(b.cost, b.index);
  }
};
