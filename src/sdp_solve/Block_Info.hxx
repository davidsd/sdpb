#pragma once

#include "Block_Cost.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <El.hpp>
#include <algorithm>
#include <filesystem>

struct MPI_Comm_Wrapper
{
  El::mpi::Comm value;
  MPI_Comm_Wrapper() = default;
  MPI_Comm_Wrapper(const MPI_Comm_Wrapper &) = delete;
  void operator=(const MPI_Comm_Wrapper &) = delete;
  ~MPI_Comm_Wrapper()
  {
    if(value != El::mpi::COMM_WORLD)
      {
        El::mpi::Free(value);
      }
  }
};

struct MPI_Group_Wrapper
{
  El::mpi::Group value;
  MPI_Group_Wrapper() = default;
  MPI_Group_Wrapper(const MPI_Group_Wrapper &) = delete;
  void operator=(const MPI_Group_Wrapper &) = delete;
  ~MPI_Group_Wrapper()
  {
    if(value != El::mpi::GROUP_NULL)
      {
        El::mpi::Free(value);
      }
  }
};

class Block_Info
{
public:
  // TODO: The filename should not be in this object.
  std::filesystem::path block_timings_filename;
  std::vector<size_t> dimensions;
  std::vector<size_t> num_points;

  std::vector<size_t> block_indices;
  MPI_Group_Wrapper mpi_group;
  MPI_Comm_Wrapper mpi_comm;

  Block_Info() = delete;
  Block_Info(const std::filesystem::path &sdp_path,
             const std::filesystem::path &checkpoint_in,
             const size_t &procs_per_node, const size_t &proc_granularity,
             const Verbosity &verbosity);
  Block_Info(const std::filesystem::path &sdp_path,
             const El::Matrix<int32_t> &block_timings,
             const size_t &procs_per_node, const size_t &proc_granularity,
             const Verbosity &verbosity);
  Block_Info(const std::vector<size_t> &matrix_dimensions,
             const size_t &procs_per_node, const size_t &proc_granularity,
             const Verbosity &verbosity);
  Block_Info(const std::vector<size_t> &matrix_dimensions,
             const Verbosity &verbosity)
      : Block_Info(matrix_dimensions, 1, 1, verbosity)
  {}
  void read_block_info(const std::filesystem::path &sdp_path);
  std::vector<Block_Cost>
  read_block_costs(const std::filesystem::path &sdp_path,
                   const std::filesystem::path &checkpoint_in);
  void
  allocate_blocks(const std::vector<Block_Cost> &block_costs,
                  const size_t &procs_per_node, const size_t &proc_granularity,
                  const Verbosity &verbosity);

  std::vector<size_t> schur_block_sizes() const
  {
    std::vector<size_t> result(num_points.size());
    for(size_t index(0); index != num_points.size(); ++index)
      {
        result[index] = num_points[index] * dimensions[index]
                        * (dimensions[index] + 1) / 2;
      }
    return result;
  }
  std::vector<size_t> bilinear_pairing_block_sizes() const
  {
    std::vector<size_t> result(2 * num_points.size());
    for(size_t index(0); index != num_points.size(); ++index)
      {
        result[2 * index] = num_points[index] * dimensions[index];
        result[2 * index + 1] = result[2 * index];
      }
    return result;
  }
  std::vector<size_t> psd_matrix_block_sizes() const
  {
    std::vector<size_t> result(2 * num_points.size());
    for(size_t index(0); index != num_points.size(); ++index)
      {
        // Need to round down (num_points+1)/2 before multiplying by
        // dim, since dim could be 2.
        result[2 * index] = dimensions[index] * ((num_points[index] + 1) / 2);
        result[2 * index + 1]
          = dimensions[index] * num_points[index] - result[2 * index];
      }
    return result;
  }
};

namespace std
{
  inline void swap(MPI_Comm_Wrapper &a, MPI_Comm_Wrapper &b)
  {
    swap(a.value, b.value);
  }
  inline void swap(MPI_Group_Wrapper &a, MPI_Group_Wrapper &b)
  {
    swap(a.value, b.value);
  }

  inline void swap(Block_Info &a, Block_Info &b)
  {
    swap(a.block_timings_filename, b.block_timings_filename);
    swap(a.dimensions, b.dimensions);
    swap(a.num_points, b.num_points);
    swap(a.block_indices, b.block_indices);
    swap(a.mpi_group, b.mpi_group);
    swap(a.mpi_comm, b.mpi_comm);
  }
}
