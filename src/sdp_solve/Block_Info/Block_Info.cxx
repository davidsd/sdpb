#include "../Block_Info.hxx"
#include "block_costs.hxx"
#include "matrix_sizes.hxx"
#include "read_block_info.hxx"
#include "sdpb_util/block_mapping/allocate_block_mapping.hxx"

namespace fs = std::filesystem;

// Constructors

Block_Info::Block_Info(const Environment &env,
                       const std::vector<size_t> &dimensions,
                       const std::vector<size_t> &num_points,
                       const std::vector<Block_Cost> &block_costs,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
    : dimensions(dimensions), num_points(num_points)
{
  ASSERT_EQUAL(dimensions.size(), num_points.size());
  ASSERT_EQUAL(dimensions.size(), block_costs.size());

  const auto print_block
    = [this](std::ostream &os, const size_t block_index) -> auto & {
    return os << block_index << "(" << this->dimensions[block_index] << ","
              << this->num_points[block_index] << ")";
  };
  allocate_block_mapping(env, block_costs, proc_granularity, print_block,
                         verbosity, mpi_group.value, mpi_comm.value,
                         block_indices);
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
  swap(a.dimensions, b.dimensions);
  swap(a.num_points, b.num_points);
  swap(a.block_indices, b.block_indices);
  swap(a.mpi_group, b.mpi_group);
  swap(a.mpi_comm, b.mpi_comm);
}
