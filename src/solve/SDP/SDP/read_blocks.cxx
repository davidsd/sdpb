#include "read_vector.hxx"
#include "../../SDP.hxx"

void read_blocks(const boost::filesystem::path &sdp_directory, SDP &sdp)
{
  boost::filesystem::ifstream block_stream(sdp_directory / "blocks");
  if(!block_stream.good())
    {
      throw std::runtime_error("Could not open '"
                               + (sdp_directory / "blocks").string() + "'");
    }
  read_vector(block_stream, sdp.dimensions);
  read_vector(block_stream, sdp.degrees);
  read_vector(block_stream, sdp.schur_block_sizes);
  read_vector(block_stream, sdp.psd_matrix_block_sizes);
  read_vector(block_stream, sdp.bilinear_pairing_block_sizes);
}
