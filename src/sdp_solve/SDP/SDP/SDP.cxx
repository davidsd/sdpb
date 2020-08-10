#include "../../SDP.hxx"

#include <boost/filesystem.hpp>

void read_blocks(const boost::filesystem::path &sdp_directory, SDP &sdp);
void compute_block_grid_mapping(const size_t &num_blocks,
                                std::vector<size_t> &block_indices);
void read_objectives(const boost::filesystem::path &sdp_directory,
                     const El::Grid &grid, El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b);
void read_bilinear_bases(
  const boost::filesystem::path &sdp_directory, const Block_Info &block_info,
  const El::Grid &grid,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases_dist);

void read_primal_objective_c(const boost::filesystem::path &sdp_directory,
                             const std::vector<size_t> &block_indices,
                             const El::Grid &grid,
                             Block_Vector &primal_objective_c);
void read_free_var_matrix(const boost::filesystem::path &sdp_directory,
                          const std::vector<size_t> &block_indices,
                          const El::Grid &grid, Block_Matrix &free_var_matrix);

SDP::SDP(const boost::filesystem::path &sdp_directory,
         const Block_Info &block_info, const El::Grid &grid)
{
  read_objectives(sdp_directory, grid, objective_const, dual_objective_b);
  read_bilinear_bases(sdp_directory, block_info, grid,
                      bilinear_bases_local, bilinear_bases_dist);
  read_primal_objective_c(sdp_directory, block_info.block_indices, grid,
                          primal_objective_c);
  read_free_var_matrix(sdp_directory, block_info.block_indices, grid,
                       free_var_matrix);
}
