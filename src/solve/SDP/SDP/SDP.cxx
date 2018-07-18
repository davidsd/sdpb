#include "read_vector.hxx"
#include "../../SDP.hxx"

#include <boost/filesystem.hpp>

SDP::SDP(const boost::filesystem::path &sdp_directory)
{
  boost::filesystem::ifstream block_stream(sdp_directory / "blocks");
  read_vector(block_stream, dimensions);
  read_vector(block_stream, degrees);
  read_vector(block_stream, schur_block_sizes);
  read_vector(block_stream, psd_matrix_block_sizes);
  read_vector(block_stream, bilinear_pairing_block_sizes);
  exit(0);
}

