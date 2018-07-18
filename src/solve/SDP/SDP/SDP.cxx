#include "../../SDP.hxx"

#include <boost/filesystem.hpp>

void read_blocks(const boost::filesystem::path &sdp_directory, SDP &sdp);
void compute_block_grid_mapping(const size_t &num_blocks,
                                std::vector<size_t> &block_indices);
void read_objectives(const boost::filesystem::path &sdp_directory,
                     El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b);
void read_bilinear_bases(
  const boost::filesystem::path &sdp_directory, const El::Grid &grid,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases_dist);

SDP::SDP(const boost::filesystem::path &sdp_directory)
    : grid(El::mpi::COMM_SELF)

{
  read_blocks(sdp_directory, *this);
  compute_block_grid_mapping(dimensions.size(), block_indices);
  read_objectives(sdp_directory, objective_const, dual_objective_b);
  read_bilinear_bases(sdp_directory, grid, bilinear_bases_local,
                      bilinear_bases_dist);
  // read_primal_objective(output_dir, primal_objective);
  // read_free_var_matrix(output_dir, free_var_matrix);

  El::mpi::Finalize();
  exit(0);
}
