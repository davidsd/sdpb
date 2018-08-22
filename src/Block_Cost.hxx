#pragma once

#include <vector>

struct Block_Cost
{
  size_t cost, index;
  Block_Cost(const size_t &Cost, const size_t &Index)
      : cost(Cost), index(Index)
  {}
  bool operator<(const Block_Cost &b) const { return cost < b.cost; }
};
