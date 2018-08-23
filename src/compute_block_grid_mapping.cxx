// Allocate blocks to MPI processes.  Blocks can be split among
// multiple processes, and multiple blocks can be assigned to a single
// process.  The algorithm ensures that all of the group's of processes
// do not cross node boundaries.

#include "Block_Cost.hxx"
#include "Block_Map.hxx"

#include <numeric>
#include <algorithm>
#include <stdexcept>

std::vector<std::vector<Block_Map>>
compute_block_grid_mapping(const size_t &procs_per_node,
                           const size_t &num_nodes,
                           const std::vector<Block_Cost> &block_costs)
{
  const double total_cost(
    std::accumulate(block_costs.begin(), block_costs.end(), 0.0,
                    [](const double &cost, const Block_Cost &element) {
                      return cost + element.cost;
                    }));
  const size_t num_procs(procs_per_node * num_nodes);
  const double average_cost(total_cost / num_procs);

  std::vector<size_t> available_procs(num_nodes, procs_per_node);

  std::vector<std::vector<Block_Map>> result;
  std::vector<Block_Map> block_map_vector;
  Block_Map block_map;
  size_t node(0);
  for(auto &block : block_costs)
    {
      if(block_map.num_procs != 0
         && block_map.num_procs * average_cost <= block_map.cost)
        {
          block_map_vector.push_back(block_map);
          block_map.clear();
        }
      if(block_map.num_procs == 0)
        {
          size_t num_procs(
            std::min(procs_per_node,
                     std::max(size_t(1), size_t(block.cost / average_cost))));
          if(num_procs > available_procs[node])
            {
              ++node;
              if(node == num_nodes)
                {
                  throw std::runtime_error(
                    "INTERNAL ERROR: Too many processors used "
                    "in the initial distribution of blocks.  Maybe try "
                    "increasing procs_per_node?");
                }
              result.push_back(block_map_vector);
              block_map_vector.clear();
            }
          block_map.num_procs = num_procs;
          available_procs[node] -= num_procs;
        }
      block_map.cost += block.cost;
      block_map.block_indices.push_back(block.index);
    }
  block_map_vector.push_back(block_map);
  result.push_back(block_map_vector);

  for(node = 0; node < num_nodes; ++node)
    {
      while(available_procs[node] != 0)
        {
          auto max_element(
            std::max_element(result[node].begin(), result[node].end()));
          if(max_element->block_indices.size() != 1)
            {
              size_t half_cost(max_element->cost / 2);
              std::vector<size_t> old_indices;
              std::swap(max_element->block_indices, old_indices);
              max_element->cost = 0;

              Block_Map new_element;
              new_element.num_procs = 1;
              --(available_procs[node]);

              for(auto &index : old_indices)
                {
                  auto cost_index(std::find_if(
                    block_costs.begin(), block_costs.end(),
                    [&](const Block_Cost &a) { return a.index == index; }));
                  if(cost_index == block_costs.end())
                    {
                      throw std::runtime_error(
                        "INTERNAL ERROR: Could not find index "
                        + std::to_string(index) + " in the block_costs array");
                    }
                  const size_t block_cost(cost_index->cost);
                  if(max_element->cost + block_cost <= half_cost)
                    {
                      max_element->block_indices.push_back(index);
                      max_element->cost += block_cost;
                    }
                  else
                    {
                      new_element.block_indices.push_back(index);
                      new_element.cost += block_cost;
                    }
                }
              result[node].push_back(new_element);
            }
          else
            {
              ++(max_element->num_procs);
              --(available_procs[node]);
            }
        }
    }
  return result;
}
