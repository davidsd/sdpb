#include "read_vector.hxx"
#include "../../SDP.hxx"

#include <boost/filesystem.hpp>

void compute_block_grid_mapping(const size_t &num_blocks,
                                std::vector<size_t> &block_indices);

SDP::SDP(const boost::filesystem::path &sdp_directory)
    : grid(El::mpi::COMM_SELF)

{
  boost::filesystem::ifstream block_stream(sdp_directory / "blocks");
  read_vector(block_stream, dimensions);
  read_vector(block_stream, degrees);
  read_vector(block_stream, schur_block_sizes);
  read_vector(block_stream, psd_matrix_block_sizes);
  read_vector(block_stream, bilinear_pairing_block_sizes);

  compute_block_grid_mapping(dimensions.size(), block_indices);
  exit(0);
}
