#include "assign_bilinear_bases_dist.hxx"
#include "set_dual_objective_b.hxx"
#include "../../SDP.hxx"

#include <boost/filesystem.hpp>

void read_blocks(const boost::filesystem::path &sdp_directory, SDP &sdp);
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
  read_bilinear_bases(sdp_directory, block_info, grid, bilinear_bases_local,
                      bilinear_bases_dist);
  read_primal_objective_c(sdp_directory, block_info.block_indices, grid,
                          primal_objective_c);
  read_free_var_matrix(sdp_directory, block_info.block_indices, grid,
                       free_var_matrix);
}

SDP::SDP(const El::BigFloat &objective_const_input,
         const std::vector<El::BigFloat> &dual_objective_b_input,
         const std::vector<std::vector<El::BigFloat>> &primal_objective_c_input,
         const std::vector<El::Matrix<El::BigFloat>> &free_var_input,
         const Block_Info &block_info, const El::Grid &grid):
  objective_const(objective_const_input)
{
  set_dual_objective_b(dual_objective_b_input, grid, dual_objective_b);

  auto &block_indices(block_info.block_indices);
  bilinear_bases_local.resize(2 * block_indices.size());

  for(size_t block(0); block != block_indices.size(); ++block)
    {
      bilinear_bases_local[2 * block].Resize(1, 1);
      bilinear_bases_local[2 * block + 1].Resize(0, 1);
      bilinear_bases_local[2 * block](0, 0) = 1;
    }
  assign_bilinear_bases_dist(bilinear_bases_local, grid, bilinear_bases_dist);

  for(size_t block(0); block != block_indices.size(); ++block)
    {
      primal_objective_c.blocks.emplace_back(
        primal_objective_c_input.at(block_indices.at(block)).size(), 1, grid);
      // TODO: copied from sdp_solve/SDP/SDP/read_primal_objective_c.cxx
      auto &primal_block(primal_objective_c.blocks.back());
      size_t local_height(primal_block.LocalHeight());
      if(primal_block.GlobalCol(0) == 0)
        {
          for(size_t row = 0; row < local_height; ++row)
            {
              size_t global_row(primal_block.GlobalRow(row));
              primal_block.SetLocal(
                row, 0,
                primal_objective_c_input.at(block_indices.at(block))
                  .at(global_row));
            }
        }
    }

  free_var_matrix.blocks.reserve(block_indices.size());
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      free_var_matrix.blocks.emplace_back(
        free_var_input.at(block_indices.at(block)).Height(),
        free_var_input.at(block_indices.at(block)).Width(), grid);
      auto &free_var_block(free_var_matrix.blocks.back());
      for(int64_t row(0); row != free_var_block.Height(); ++row)
        for(int64_t column(0); column != free_var_block.Width(); ++column)
          {
            if(free_var_block.IsLocal(row, column))
              {
                free_var_block.SetLocal(
                  free_var_block.LocalRow(row),
                  free_var_block.LocalCol(column),
                  free_var_input.at(block_indices.at(block))(row, column));
              }
          }
    }
}
