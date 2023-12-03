// Allocate blocks to MPI processes.  Blocks can be split among
// multiple processes, and multiple blocks can be assigned to a single
// process.  The algorithm ensures that all of the group's of processes
// do not cross node boundaries.

// With multiple nodes, this turns into a bin-packing problem.  The
// major difference from the classical bin-packing problem is that we
// have a fixed set of bins.  The capacity of each node is also
// somewhat elastic, in that we can cram more blocks into a node if
// they will not fit anywhere else.

// This algorithm starts by assigning blocks to nodes, where the
// number of procs for a given block is
//
//   num_procs=min(block_cost/average_cost,1)
//
// rounded down.  The block is assigned to the node with the most
// remaining spots.  This is a Worst Fit First method, which is
// usually not the first choice for bin-packing.  We do this because
// this method usually overloads these multi-cpu blocks.  Putting them
// on separate nodes means that, if there are any left-over procs,
// they can be separately allocated to these high cost blocks.
//
// It can happen that it will not fit on any node, because each node
// only has a few procs left.  In that case, the block will get
// assigned to whatever node has the most leftover procs.
//
// This will usually overload some procs, but guarantees that all
// blocks will get allocated somewhere.  This usually results in some
// left over procs on nodes.  We assign these left over procs to the
// highest cost block_map, where the cost is the total cost of all
// blocks divided by by the number of procs assigned to that
// block_map.
//
// If the block_map has multiple blocks, then instead of assiging more
// procs, we create a new block_map.  The blocks are then split
// between the old and new block_maps.
//
// This algorithm is certainly far from optimal.
//
// 1) When block_maps are split, they are not rebalanced amongst all
// the block_maps.
//
// 2) When large blocks are forced to fit into a node, there is no
// sharing of procs between the existing block_maps and the new entry.

#include "sdp_solve/Block_Cost.hxx"
#include "Block_Map.hxx"

#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <limits>

std::vector<std::vector<Block_Map>>
compute_block_grid_mapping(const size_t &procs_per_node,
                           const size_t &num_nodes,
                           const std::vector<Block_Cost> &block_costs)
{
  // We do computations in integers to make sure that the results are
  // the same on different processers.
  const size_t total_cost(
    std::accumulate(block_costs.begin(), block_costs.end(), size_t(0),
                    [](const size_t &cost, const Block_Cost &element) {
                      return cost + element.cost;
                    }));
  const size_t num_procs(procs_per_node * num_nodes);
  std::vector<size_t> available_procs(num_nodes, procs_per_node);

  std::vector<std::vector<Block_Map>> result(num_nodes);

  auto multi_proc_end(std::find_if(
    block_costs.begin(), block_costs.end(), [&](const Block_Cost &block) {
      return num_procs * block.cost <= total_cost;
    }));

  size_t remaining_cost(total_cost), remaining_procs(num_procs);
  for(auto block(block_costs.begin()); block != multi_proc_end; ++block)
    {
      // Always add block_map's to the node with the most available
      // procs.  This is Worst Fit First.

      size_t max_available_node(std::distance(
        available_procs.begin(),
        std::max_element(available_procs.begin(), available_procs.end())));

      size_t procs_for_block = std::min(
        available_procs[max_available_node],
        std::max(size_t(1), size_t(block->cost * num_procs / total_cost)));
      Block_Map block_map(procs_for_block, block->cost, {block->index});
      result[max_available_node].push_back(block_map);
      available_procs[max_available_node] -= procs_for_block;
      remaining_cost -= block->cost;
      remaining_procs -= procs_for_block;
    }

  if(remaining_procs != num_procs)
    {
      size_t required_procs((remaining_cost * num_procs + (total_cost - 1))
                            / total_cost);
      size_t extra_procs(remaining_procs - required_procs);
      while(extra_procs > 0)
        {
          std::vector<Block_Map>::iterator max_element(result.front().end());
          size_t max_node;
          for(size_t node = 0; node < result.size(); ++node)
            {
              if(available_procs[node] > 0 && !result[node].empty())
                {
                  auto max_element_on_node(std::max_element(
                    result[node].begin(), result[node].end()));
                  if(max_element == result.front().end()
                     || *max_element < *max_element_on_node)
                    {
                      max_element = max_element_on_node;
                      max_node = node;
                    }
                }
            }
          if(max_element != result.front().end())
            {
              max_element->num_procs += 1;
              --available_procs[max_node];
              --extra_procs;
            }
          else
            {
              break;
            }
        }
    }
  std::vector<std::vector<Block_Map>> available_block_maps(num_nodes);
  for(size_t node = 0; node < num_nodes; ++node)
    {
      available_block_maps[node].resize(available_procs[node]);
      for(auto &block : available_block_maps[node])
        {
          block.num_procs = 1;
        }
    }

  for(auto block(multi_proc_end); block != block_costs.end(); ++block)
    {
      size_t min_cost(std::numeric_limits<size_t>::max());
      auto min_block(available_block_maps.at(0).end());
      for(size_t node(0); node < num_nodes; ++node)
        {
          if(!available_block_maps[node].empty())
            {
              auto block = std::min_element(available_block_maps[node].begin(),
                                            available_block_maps[node].end());
              if(block->cost < min_cost)
                {
                  min_block = block;
                  min_cost = block->cost;
                }
            }
        }
      if(min_block==available_block_maps.at(0).end())
        {
          throw std::runtime_error("INTERNAL ERROR: Unable to find any "
                                   "free processors for remaining blocks");
        }
      min_block->cost += block->cost;
      min_block->block_indices.push_back(block->index);
    }
  for(size_t node = 0; node < num_nodes; ++node)
    for(auto &available : available_block_maps[node])
      {
        result[node].push_back(available);
      }
  return result;
}
