#include "Block_Cost.hxx"
#include "Block_Map.hxx"

#include <numeric>
#include <algorithm>
#include <stdexcept>

std::vector<Block_Map>
compute_block_grid_mapping(const size_t &num_procs,
                           const std::vector<Block_Cost> &block_costs)
{
  double total_cost(
    std::accumulate(block_costs.begin(), block_costs.end(), 0.0,
                    [](const double &cost, const Block_Cost &element) {
                      return cost + element.cost;
                    }));
  double average_cost(total_cost / num_procs);

  size_t available_procs(num_procs);
  Block_Map block_map;

  std::vector<Block_Map> result;
  for(auto &block: block_costs)
    {
      if(block_map.num_procs != 0
         && block_map.num_procs * average_cost <= block_map.cost)
        {
          result.push_back(block_map);
          block_map.clear();
        }
      if(block_map.num_procs == 0)
        {
          size_t num_procs(
            std::max(1.0, block.cost / average_cost));
          if(num_procs > available_procs)
            {
              throw std::runtime_error(
                "INTERNAL ERROR: Too many processors used "
                "in the initial distribution of blocks");
            }
          block_map.num_procs = num_procs;
          available_procs -= num_procs;
        }
      block_map.cost += block.cost;
      block_map.block_indices.push_back(block.index);
    }
  result.push_back(block_map);

  while(available_procs != 0)
    {
      std::sort(result.begin(), result.end());
      ++result.front().num_procs;
      --available_procs;
    }
  return result;
}
