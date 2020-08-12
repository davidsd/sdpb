#pragma once

#include "read_vector.hxx"
#include "Verbosity.hxx"
#include "../Block_Cost.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>

#include <algorithm>

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
  boost::filesystem::path block_timings_filename;
  size_t file_num_procs;
  std::vector<std::vector<size_t>> file_block_indices;

  // dimensions[j] = m_j  (0 <= j < J)
  std::vector<size_t> dimensions;

  // degrees[j] = d_j  (0 <= j < J)
  std::vector<size_t> degrees;

  std::vector<size_t> schur_block_sizes,
    // Dimensions of the blocks of X,Y
    // psd_matrix_block_sizes[b] = (delta_b+1)*m_j = length(v_{b,*})*m_j
    // (0 <= b < bMax)
    psd_matrix_block_sizes,
    // Dimensions of the bilinear pairing matrices U and V
    // bilinear_pairing_block_sizes[b] = (d_j + 1)*m_j
    // (0 <= b < bMax)
    bilinear_pairing_block_sizes;

  std::vector<size_t> block_indices;
  MPI_Group_Wrapper mpi_group;
  MPI_Comm_Wrapper mpi_comm;

  Block_Info() = delete;
  Block_Info(const boost::filesystem::path &sdp_directory,
             const boost::filesystem::path &checkpoint_in,
             const size_t &procs_per_node, const size_t &proc_granularity,
             const Verbosity &verbosity);
  Block_Info(const boost::filesystem::path &sdp_directory,
             const El::Matrix<int32_t> &block_timings,
             const size_t &procs_per_node, const size_t &proc_granularity,
             const Verbosity &verbosity);
  Block_Info(const std::vector<size_t> &matrix_dimensions,
             const size_t &procs_per_node, const size_t &proc_granularity,
             const Verbosity &verbosity);
  void read_block_info(const boost::filesystem::path &sdp_directory);
  std::vector<Block_Cost>
  read_block_costs(const boost::filesystem::path &sdp_directory,
                   const boost::filesystem::path &checkpoint_in);
  void
  allocate_blocks(const std::vector<Block_Cost> &block_costs,
                  const size_t &procs_per_node, const size_t &proc_granularity,
                  const Verbosity &verbosity);
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
    swap(a.file_num_procs, b.file_num_procs);
    swap(a.file_block_indices, b.file_block_indices);
    swap(a.dimensions, b.dimensions);
    swap(a.degrees, b.degrees);
    swap(a.schur_block_sizes, b.schur_block_sizes);
    swap(a.psd_matrix_block_sizes, b.psd_matrix_block_sizes);
    swap(a.bilinear_pairing_block_sizes, b.bilinear_pairing_block_sizes);
    swap(a.block_indices, b.block_indices);
    swap(a.mpi_group, b.mpi_group);
    swap(a.mpi_comm, b.mpi_comm);
  }
}
