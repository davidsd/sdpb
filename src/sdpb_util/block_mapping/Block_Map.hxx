#pragma once

#include <limits>
#include <vector>

struct Block_Map
{
  // Either of num_procs>1 or block_indices.size()>1 can be true, but
  // not both.
  size_t num_procs = 0;
  size_t cost = 0;
  std::vector<size_t> block_indices;

  Block_Map() = default;
  Block_Map(const size_t &Num_procs, const size_t &Cost,
            const std::initializer_list<size_t> &Block_indices)
      : num_procs(Num_procs), cost(Cost), block_indices(Block_indices)
  {}
  void clear()
  {
    num_procs = cost = 0;
    block_indices.clear();
  }
  // Sort by average cost
  bool operator<(const Block_Map &b) const
  {
    return std::make_tuple(cost * b.num_procs, first_index())
           < std::make_tuple(b.cost * num_procs, b.first_index());
  }

private:
  [[nodiscard]] size_t first_index() const
  {
    return block_indices.empty() ? std::numeric_limits<size_t>::max()
                                 : block_indices.at(0);
  }
};
