#pragma once

#include "sdpb_util/assert.hxx"
#include "sdpb_util/block_mapping/Block_Cost.hxx"

#include <vector>

inline std::vector<Block_Cost>
read_block_costs_from_timings(const std::filesystem::path &block_timings_path,
                              const size_t num_blocks)
{
  ASSERT(exists(block_timings_path));
  std::vector<Block_Cost> result;
  result.reserve(num_blocks);
  size_t index(0), cost;
  std::ifstream costs(block_timings_path);
  costs >> cost;
  while(costs.good())
    {
      result.emplace_back(cost, index);
      ++index;
      costs >> cost;
    }
  if(result.size() != num_blocks)
    {
      RUNTIME_ERROR("Incompatible number of entries in ", block_timings_path,
                    ": expected ", num_blocks, " but found ", result.size());
    }
  return result;
}
