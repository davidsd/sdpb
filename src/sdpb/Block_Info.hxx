#pragma once

#include "read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>

class Block_Info
{
public:
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
  El::mpi::Group mpi_group;
  El::mpi::Comm mpi_comm;

  Block_Info() = delete;
  Block_Info(const boost::filesystem::path &sdp_directory,
             const boost::filesystem::path &checkpoint_in,
             const size_t &procs_per_node);
};
