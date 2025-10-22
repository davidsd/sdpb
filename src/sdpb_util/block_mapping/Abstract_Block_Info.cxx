#include "Abstract_Block_Info.hxx"
#include "allocate_block_mapping.hxx"
#include "sdpb_util/assert.hxx"

namespace fs = std::filesystem;

// Constructors

Abstract_Block_Info::Abstract_Block_Info(
  const Environment &env, const std::vector<Block_Cost> &block_costs,
  const size_t &proc_granularity,
  const std::function<std::ostream &(std::ostream &, size_t)> &print_block,
  const Verbosity &verbosity)
    : block_mapping(allocate_block_mapping(
        env, block_costs, proc_granularity, print_block, verbosity,
        mpi_group.value, mpi_comm.value, block_indices)),
      node_comm(env.comm_shared_mem),
      node_index(env.node_index())
{
  // Some sanity checks
  const int node_rank = node_comm.Rank();
  const int num_nodes = env.num_nodes();
  const int num_node_ranks = node_comm.Size();
  ASSERT_EQUAL(num_nodes, block_mapping.mapping.size());
  ASSERT_EQUAL(num_nodes, block_mapping.node_rank_to_group.size());
  ASSERT_EQUAL(num_node_ranks,
               block_mapping.node_rank_to_group.at(node_index).size());
  const auto group_index = node_group_index();
  ASSERT(
    this->block_indices
      == block_mapping.mapping.at(node_index).at(group_index).block_indices,
    "Wrong block indices", DEBUG_STRING(node_index), DEBUG_STRING(node_rank),
    DEBUG_STRING(group_index));
}

size_t Abstract_Block_Info::node_group_index() const
{
  return block_mapping.node_rank_to_group.at(node_index).at(node_comm.Rank());
}

void swap(Abstract_Block_Info &a, Abstract_Block_Info &b) noexcept
{
  using std::swap;
  swap(a.block_indices, b.block_indices);
  swap(a.mpi_group, b.mpi_group);
  swap(a.mpi_comm, b.mpi_comm);
  swap(a.block_mapping, b.block_mapping);
  swap(a.node_comm, b.node_comm);
  swap(a.node_index, b.node_index);
}
