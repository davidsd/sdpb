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
size_t Abstract_Block_Info::num_nodes() const
{
  return block_mapping.mapping.size();
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

void Abstract_Block_Info::allocate_single_rank_blocks(const std::vector<int> &owner)
{
  const int rank = El::mpi::Rank(), num_ranks = El::mpi::Size();
  ASSERT_EQUAL(owner.size(), block_mapping.block_locations.size());

  // Each rank's place (node, rank on the node), from every rank, so that all
  // ranks build the same mapping whatever the launcher's placement, and every
  // check below sees the same data on every rank
  std::vector<int> place(2 * size_t(num_ranks));
  const int mine[2] = {int(node_index), node_comm.Rank()};
  MPI_Allgather(mine, 2, MPI_INT, place.data(), 2, MPI_INT, MPI_COMM_WORLD);

  // Node -> group (= rank on the node) -> its blocks
  const size_t num_nodes = this->num_nodes();
  std::vector<std::vector<Block_Map>> mapping(num_nodes);
  auto group_of = [&](const int r) -> Block_Map & {
    return mapping.at(size_t(place[2 * r])).at(size_t(place[2 * r + 1]));
  };
  for(int r = 0; r < num_ranks; ++r)
    {
      const size_t node = size_t(place[2 * r]),
                   node_rank = size_t(place[2 * r + 1]);
      ASSERT(node < num_nodes, "rank ", r, " is on node ", node, " of ",
             num_nodes);
      if(mapping[node].size() <= node_rank)
        mapping[node].resize(node_rank + 1);
      mapping[node][node_rank].num_procs = 1;
    }
  for(size_t block = 0; block < owner.size(); ++block)
    {
      const int r = owner[block];
      ASSERT(r >= 0 && r < num_ranks, "Invalid GPU block owner: ", r,
             DEBUG_STRING(block));
      group_of(r).block_indices.push_back(block);
    }
  for(int r = 0; r < num_ranks; ++r)
    ASSERT(!group_of(r).block_indices.empty(), "GPU rank ", r,
           " has no chunks; lower SDPB_GPU_CHUNK_BLOCKS");

  MPI_Comm_Wrapper comm;
  El::mpi::Split(El::mpi::COMM_WORLD, rank, 0, comm.value);
  MPI_Group_Wrapper group;
  El::mpi::CommGroup(comm.value, group.value);
  // The old communicator and group are freed with comm and group
  swap(mpi_comm, comm);
  swap(mpi_group, group);
  block_mapping = Block_Mapping(mapping);
  block_indices = block_mapping.mapping.at(node_index)
                    .at(node_group_index())
                    .block_indices;
}
