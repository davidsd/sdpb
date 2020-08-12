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

SDP::SDP(const std::vector<El::BigFloat> &objectives,
         const std::vector<El::BigFloat> &normalization,
         const std::vector<El::BigFloat> &prefactor,
         const std::vector<std::vector<El::BigFloat>> &primal_objective_c_input,
         const std::vector<El::Matrix<El::BigFloat>> &free_var_input,
         const Block_Info &block_info, const El::Grid &grid)
{
  // TODO: This is duplicated from sdp2input/write_output/write_output.cxx
  auto max_normalization(normalization.begin());
  for(auto n(normalization.begin()); n != normalization.end(); ++n)
    {
      if(Abs(*n) > Abs(*max_normalization))
        {
          max_normalization = n;
        }
    }
  size_t max_index(std::distance(normalization.begin(), max_normalization));

  objective_const = objectives.at(max_index) / normalization.at(max_index);

  std::vector<El::BigFloat> offset_dual_objective_b;
  offset_dual_objective_b.reserve(normalization.size() - 1);
  for(size_t index = 0; index < normalization.size(); ++index)
    {
      if(index != max_index)
        {
          offset_dual_objective_b.push_back(
            objectives.at(index) - normalization.at(index) * objective_const);
        }
    }

  set_dual_objective_b(offset_dual_objective_b, grid, dual_objective_b);

  auto &block_indices(block_info.block_indices);
  bilinear_bases_local.resize(2 * block_indices.size());

  for(size_t block(0); block != block_indices.size(); ++block)
    {
      bilinear_bases_local[2 * block].Resize(1, 1);
      bilinear_bases_local[2 * block](0, 0) = Sqrt(prefactor.at(block));
    }
  assign_bilinear_bases_dist(bilinear_bases_local, grid, bilinear_bases_dist);

  for(size_t index(0); index != block_indices.size(); ++index)
    {
      primal_objective_c.blocks.emplace_back(
        primal_objective_c_input.at(index).size(), 1, grid);
      // TODO: copied from sdp_solve/SDP/SDP/read_primal_objective_c.cxx
      auto &block(primal_objective_c.blocks.back());
      size_t local_height(block.LocalHeight());
      if(block.GlobalCol(0) == 0)
        {
          for(size_t row = 0; row < local_height; ++row)
            {
              size_t global_row(block.GlobalRow(row));
              block.SetLocal(
                row, 0, primal_objective_c_input.at(index).at(global_row));
            }
        }
    }

  free_var_matrix.blocks.reserve(block_indices.size());
  for(size_t index(0); index != block_indices.size(); ++index)
    {
      free_var_matrix.blocks.emplace_back(free_var_input.at(index).Height(),
                                          free_var_input.at(index).Width(),
                                          grid);
      auto &block(free_var_matrix.blocks.back());
      for(int64_t row(0); row != block.Height(); ++row)
        for(int64_t column(0); column != block.Height(); ++column)
          {
            if(block.IsLocal(row, column))
              {
                block.SetLocal(block.LocalRow(row), block.LocalCol(column),
                               free_var_input.at(index)(row, column));
              }
          }
    }
}
