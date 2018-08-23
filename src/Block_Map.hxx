#pragma once

#include <vector>

struct Block_Map
{
  // Either of num_procs>1 or block_indices.size()>1 can be true, but
  // not both.
  size_t num_procs=0;
  size_t cost=0;
  std::vector<size_t> block_indices;
  void clear()
  {
    num_procs=cost=0;
    block_indices.clear();
  }
  // Sort by average cost
  bool operator<(const Block_Map &b) const
  {
    return cost*b.num_procs < b.cost*num_procs;
  }
};
