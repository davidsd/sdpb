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
         const Block_Info &block_info, const El::Grid &grid)
    : yp_to_y(dual_objective_b_input.size(), dual_objective_b_input.size(),
              grid),
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
      const size_t local_height(primal_block.LocalHeight());
      if(primal_block.GlobalCol(0) == 0)
        {
          for(size_t row = 0; row < local_height; ++row)
            {
              const size_t global_row(primal_block.GlobalRow(row));
              primal_block.SetLocal(
                row, 0,
                primal_objective_c_input.at(block_indices.at(block))
                  .at(global_row));
            }
        }
    }

  const int64_t B_Height([&]() {
    int64_t result(0);
    for(auto &block : primal_objective_c.blocks)
      {
        result += block.Height();
      }
    return result;
  }());
  const int64_t B_Width(dual_objective_b.Height());
  El::DistMatrix<El::BigFloat> B(B_Height, B_Width, grid);

  std::cout << "B: " << B.Height() << " " << B.Width() << "\n" << std::flush;

  int64_t row_block(0);
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      const int64_t block_height(
        free_var_input.at(block_indices.at(block)).Height()),
        block_width(free_var_input.at(block_indices.at(block)).Width());
      for(int64_t row_offset(0); row_offset != block_height; ++row_offset)
        {
          const size_t row(row_block + row_offset);
          for(int64_t column(0); column != block_width; ++column)
            {
              if(B.IsLocal(row, column))
                {
                  B.SetLocal(B.LocalRow(row), B.LocalCol(column),
                             free_var_input.at(block_indices.at(block))(
                               row_offset, column));
                }
            }
        }
      row_block += block_height;
    }

  El::DistMatrix<El::BigFloat> U(grid), s(grid), V(grid);
  El::SVD(B, U, s, V);

  for(int64_t row(0); row < s.LocalHeight(); ++row)
    {
      for(int64_t column(0); column < s.LocalWidth(); ++column)
        {
          s.SetLocal(row, column, 1 / s.GetLocal(row, column));
        }
    }
  El::Copy(V, yp_to_y);
  El::DiagonalScale(El::LeftOrRight::RIGHT, El::Orientation::NORMAL, s,
                    yp_to_y);
  
  El::DistMatrix<El::BigFloat> b_new(s.Height(), s.Width(), grid);
  El::Copy(dual_objective_b, b_new);
  El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(1.0), yp_to_y, b_new,
           El::BigFloat(0.0), dual_objective_b);

  El::Print(s, "s");
  El::Print(b_new, "\nb");
  El::Print(dual_objective_b, "\nb_prime");
  std::cout << "\n" << std::flush;

  free_var_matrix.blocks.reserve(block_indices.size());
  int64_t global_row(0);
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      free_var_matrix.blocks.emplace_back(
        free_var_input.at(block_indices.at(block)).Height(),
        free_var_input.at(block_indices.at(block)).Width(), grid);
      auto &free_var_block(free_var_matrix.blocks.back());
      for(int64_t row(0); row != free_var_block.Height(); ++row)
        {
          for(int64_t column(0); column != free_var_block.Width(); ++column)
            {
              if(free_var_block.IsLocal(row, column))
                {
                  free_var_block.SetLocal(free_var_block.LocalRow(row),
                                          free_var_block.LocalCol(column),
                                          U.Get(global_row, column));
                }
            }
          ++global_row;
        }
    }
}
