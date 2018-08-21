#pragma once

#include "read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>

class Block_Info
{
public:
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
  El::mpi::Group group;
  El::mpi::Comm comm;

  Block_Info() = delete;
  Block_Info(const boost::filesystem::path &sdp_directory)
  {
    boost::filesystem::ifstream block_stream(sdp_directory / "blocks");
    if(!block_stream.good())
      {
        throw std::runtime_error("Could not open '"
                                 + (sdp_directory / "blocks").string() + "'");
      }
    read_vector(block_stream, dimensions);
    read_vector(block_stream, degrees);
    read_vector(block_stream, schur_block_sizes);
    read_vector(block_stream, psd_matrix_block_sizes);
    read_vector(block_stream, bilinear_pairing_block_sizes);

    compute_block_grid_mapping();
  }
  void compute_block_grid_mapping();
};
