#include "../Block_Info.hxx"
#include "block_costs.hxx"
#include "matrix_sizes.hxx"
#include "read_block_info.hxx"

namespace fs = std::filesystem;

// Constructors

Block_Info::Block_Info(const Environment &env,
                       const std::vector<size_t> &dimensions,
                       const std::vector<size_t> &num_points,
                       const std::vector<Block_Cost> &block_costs,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
    : Abstract_Block_Info(
        env, block_costs, proc_granularity,
        [&dimensions, &num_points](std::ostream &os,
                                   const size_t block_index) -> auto & {
          return os << block_index << "(" << dimensions[block_index] << ","
                    << num_points[block_index] << ")";
        },
        verbosity),
      dimensions(dimensions),
      num_points(num_points)
{
  ASSERT_EQUAL(dimensions.size(), num_points.size());
  ASSERT_EQUAL(dimensions.size(), block_costs.size());
}

void Block_Info::allocate_single_rank_blocks(const std::vector<int> &owner)
{
  const int rank = El::mpi::Rank(), num_ranks = El::mpi::Size();
  ASSERT_EQUAL(owner.size(), dimensions.size());

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

Block_Info
Block_Info::create(const Environment &env, const fs::path &sdp_path,
                   const fs::path &block_timings_path,
                   const size_t &proc_granularity, const Verbosity &verbosity)
{
  std::vector<size_t> dimensions;
  std::vector<size_t> num_points;
  read_block_info(sdp_path, dimensions, num_points);
  if(block_timings_path.empty())
    {
      const auto dual_dimension = read_dual_dimension(env, sdp_path);
      return create(env, dimensions, num_points, dual_dimension,
                    proc_granularity, verbosity);
    }

  const auto block_costs
    = read_block_costs_from_timings(block_timings_path, dimensions.size());
  return Block_Info(env, dimensions, num_points, block_costs, proc_granularity,
                    verbosity);
}

Block_Info
Block_Info::create(const Environment &env, const fs::path &sdp_path,
                   const El::Matrix<int32_t> &block_timings,
                   const size_t &proc_granularity, const Verbosity &verbosity)
{
  std::vector<size_t> dimensions;
  std::vector<size_t> num_points;
  read_block_info(sdp_path, dimensions, num_points);
  std::vector<Block_Cost> block_costs;
  ASSERT_EQUAL(block_timings.Width(), 1);
  for(int64_t block = 0; block < block_timings.Height(); ++block)
    {
      block_costs.emplace_back(block_timings(block, 0), block);
    }
  return Block_Info(env, dimensions, num_points, block_costs, proc_granularity,
                    verbosity);
}
Block_Info
Block_Info::create(const Environment &env,
                   const std::vector<size_t> &dimensions,
                   const std::vector<size_t> &num_points,
                   const size_t &dual_dimension,
                   const size_t &proc_granularity, const Verbosity &verbosity)
{
  const auto block_costs
    = block_costs_from_dimensions(dimensions, num_points, dual_dimension);
  return Block_Info(env, dimensions, num_points, block_costs, proc_granularity,
                    verbosity);
}

// Matrix sizes

size_t Block_Info::get_schur_block_size(const size_t index) const
{
  return ::get_schur_block_size(dimensions, num_points, index);
}
std::vector<size_t> Block_Info::schur_block_sizes() const
{
  return ::schur_block_sizes(dimensions, num_points);
}
size_t Block_Info::get_bilinear_pairing_block_size(const size_t index,
                                                   const size_t parity) const
{
  return ::get_bilinear_pairing_block_size(dimensions, num_points, index,
                                           parity);
}
std::vector<size_t> Block_Info::bilinear_pairing_block_sizes() const
{
  return ::bilinear_pairing_block_sizes(dimensions, num_points);
}
size_t Block_Info::get_psd_matrix_block_size(const size_t index,
                                             const size_t parity) const
{
  return ::get_psd_matrix_block_size(dimensions, num_points, index, parity);
}
std::vector<size_t> Block_Info::psd_matrix_block_sizes() const
{
  return ::psd_matrix_block_sizes(dimensions, num_points);
}
size_t Block_Info::get_bilinear_bases_height(const size_t index,
                                             const size_t parity) const
{
  // see Dual_Constraint_Group code:
  const size_t degree = num_points.at(index) - 1;
  return (degree + parity) / 2 + 1 - parity;
}
size_t
Block_Info::get_bilinear_bases_width(const size_t index, const size_t) const
{
  return num_points.at(index);
}
void swap(Block_Info &a, Block_Info &b) noexcept
{
  using std::swap;
  swap(static_cast<Abstract_Block_Info &>(a),
       static_cast<Abstract_Block_Info &>(b));
  swap(a.dimensions, b.dimensions);
  swap(a.num_points, b.num_points);
}
