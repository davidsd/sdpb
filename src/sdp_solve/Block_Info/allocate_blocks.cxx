#include "../Block_Info.hxx"
#include "sdpb_util/block_mapping/allocate_block_mapping.hxx"

#include <ostream>

void Block_Info::allocate_blocks(const Environment &env,
                                 const std::vector<Block_Cost> &block_costs,
                                 const size_t &proc_granularity,
                                 const Verbosity &verbosity)
{
  const auto print_block
    = [this](std::ostream &os, const size_t block_index) -> auto & {
    return os << block_index << "(" << dimensions[block_index] << ","
              << num_points[block_index] << ")";
  };
  allocate_block_mapping(env, block_costs, proc_granularity, print_block,
                         verbosity, mpi_group.value, mpi_comm.value,
                         block_indices);
}
