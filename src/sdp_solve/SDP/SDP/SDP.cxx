#include "assign_bilinear_bases_dist.hxx"
#include "../../SDP.hxx"
#include "set_bases_blocks.hxx"

#include <boost/filesystem.hpp>

#include <iterator>

void read_blocks(const boost::filesystem::path &sdp_directory,
                 const El::Grid &grid, const Block_Info &block_info, SDP &sdp);
void read_objectives(const boost::filesystem::path &sdp_directory,
                     const El::Grid &grid, El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b);

SDP::SDP(const boost::filesystem::path &sdp_directory,
         const Block_Info &block_info, const El::Grid &grid)
{
  read_objectives(sdp_directory, grid, objective_const, dual_objective_b);
  read_blocks(sdp_directory, grid, block_info, *this);
}

SDP::SDP(const El::BigFloat &objective_const_input,
         const std::vector<El::BigFloat> &dual_objective_b_input,
         const std::vector<std::vector<El::BigFloat>> &primal_objective_c_input,
         const std::vector<El::Matrix<El::BigFloat>> &free_var_input,
         const Block_Info &block_info, const El::Grid &grid)
    : dual_objective_b(dual_objective_b_input.size(), 1, grid),
      yp_to_y(dual_objective_b_input.size(), dual_objective_b_input.size(),
              grid),
      objective_const(objective_const_input)
{
  std::vector<size_t> block_offsets(primal_objective_c_input.size() + 1, 0);
  for(size_t p(0); p < primal_objective_c_input.size(); ++p)
    {
      block_offsets[p + 1]
        = block_offsets[p] + primal_objective_c_input[p].size();
    }

  auto &block_indices(block_info.block_indices);
  std::vector<El::Matrix<El::BigFloat>> bilinear_bases_local(
    2 * block_indices.size());
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      bilinear_bases_local[2 * block].Resize(1, 1);
      bilinear_bases_local[2 * block + 1].Resize(0, 1);
      bilinear_bases_local[2 * block](0, 0) = 1;
    }
  assign_bilinear_bases_dist(bilinear_bases_local, grid, bilinear_bases);
  set_bases_blocks(block_info, bilinear_bases_local, bases_blocks, grid);

  const El::BigFloat primal_c_scale([&]() {
    El::BigFloat max_c(0.0);
    for(size_t block(0); block != block_indices.size(); ++block)
      {
        for(auto &element :
            primal_objective_c_input.at(block_indices.at(block)))
          {
            max_c = std::max(max_c, El::Abs(element));
          }
      }
    max_c = El::mpi::AllReduce(max_c, El::mpi::MAX, El::mpi::COMM_WORLD);
    if(max_c == El::BigFloat(0.0))
      {
        return El::BigFloat(1.0);
      }
    return 1 / max_c;
  }());

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
                primal_c_scale
                  * primal_objective_c_input.at(block_indices.at(block))
                      .at(global_row));
            }
        }
    }

  // Setup B
  const int64_t B_Height(block_offsets.back()),
    B_Width(dual_objective_b_input.size());
  El::Grid global_grid;
  El::DistMatrix<El::BigFloat> B(B_Height, B_Width, global_grid);

  for(int64_t row(0); row != B.LocalHeight(); ++row)
    {
      const size_t global_row(B.GlobalRow(row));
      const auto upper_iterator(std::upper_bound(
        block_offsets.begin(), block_offsets.end(), global_row));
      const size_t block_index(
        std::distance(block_offsets.begin(), upper_iterator) - 1);
      const size_t block_row(global_row - block_offsets[block_index]);
      for(int64_t column(0); column != B.LocalWidth(); ++column)
        {
          const size_t global_column(B.GlobalCol(column));
          B.SetLocal(
            row, column,
            primal_c_scale
              * free_var_input.at(block_index)(block_row, global_column));
        }
    }

  // Compute SVD of B
  //   B = U s V^T
  // Define
  //   b''(j) = Sum(b(m) * V^T(j,m)/s(j), m)
  //   B'' = U
  // This gives
  //   minimize b''(j) y''(j) + c_0
  // with constraints
  //   B''(i,j) y''(j) <= c(i)
  // This implies
  //   y(m) = Sum((V^T(l,m)/s(l)) * y''(l), l)
  // and so to convert back to y
  //   y(m) = Sum(yp_to_y(m,l) * y''(l), l)
  //   yp_to_y(m,l) = V^T(l,m)/s(l)
  El::DistMatrix<El::BigFloat> U(B.Grid());
  {
    El::DistMatrix<El::BigFloat> temp(B.Grid()), V(B.Grid()),
      dual_objective_b_global(dual_objective_b_input.size(), 1, B.Grid());
    // SVD returns U, s, and V
    El::SVD(B, U, temp, V);

    El::DistMatrix<El::BigFloat> V_s(V);
    El::DiagonalSolve(El::LeftOrRight::RIGHT, El::Orientation::NORMAL, temp,
                      V_s);

    for(int64_t row(0); row < temp.LocalHeight(); ++row)
      {
        const int64_t global_row(temp.GlobalRow(row));
        for(int64_t column(0); column < temp.LocalWidth(); ++column)
          {
            temp.SetLocal(row, column, dual_objective_b_input.at(global_row));
          }
      }

    El::Zero(dual_objective_b_global);
    El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(1.0), V_s, temp,
             El::BigFloat(0.0), dual_objective_b_global);

    const El::BigFloat b_scale([&]() {
      El::BigFloat max(El::MaxAbs(dual_objective_b_global));
      if(max == El::BigFloat(0.0))
        {
          return El::BigFloat(1.0);
        }
      return 1 / max;
    }());

    El::DistMatrix<El::BigFloat, El::STAR, El::STAR> dual_objective_b_star(
      dual_objective_b_global);
    for(int64_t row(0); row < dual_objective_b.LocalHeight(); ++row)
      {
        const int64_t global_row(dual_objective_b.GlobalRow(row));
        for(int64_t column(0); column < dual_objective_b.LocalWidth();
            ++column)
          {
            const int64_t global_column(dual_objective_b.GlobalCol(column));
            dual_objective_b.SetLocal(
              row, column,
              b_scale
                * dual_objective_b_star.GetLocal(global_row, global_column));
          }
      }

    El::Fill(temp, El::BigFloat(1.0));
    El::DiagonalScale(El::LeftOrRight::LEFT, El::Orientation::NORMAL, temp,
                      V_s);
    El::DistMatrix<El::BigFloat, El::STAR, El::STAR> yp_to_y_star(V_s);
    for(int64_t row(0); row < yp_to_y.LocalHeight(); ++row)
      {
        const int64_t global_row(yp_to_y.GlobalRow(row));
        for(int64_t column(0); column < yp_to_y.LocalWidth(); ++column)
          {
            const int64_t global_column(yp_to_y.GlobalCol(column));
            yp_to_y.SetLocal(row, column,
                             yp_to_y_star.GetLocal(global_row, global_column));
          }
      }
  }
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> U_star(U);

  free_var_matrix.blocks.reserve(block_indices.size());
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      const size_t block_index(block_indices.at(block));
      free_var_matrix.blocks.emplace_back(
        free_var_input.at(block_index).Height(),
        free_var_input.at(block_index).Width(), grid);
      auto &free_var_block(free_var_matrix.blocks.back());
      for(int64_t row(0); row != free_var_block.Height(); ++row)
        {
          const int64_t global_row(block_offsets.at(block_index) + row);
          for(int64_t column(0); column != free_var_block.Width(); ++column)
            {
              if(free_var_block.IsLocal(row, column))
                {
                  free_var_block.SetLocal(free_var_block.LocalRow(row),
                                          free_var_block.LocalCol(column),
                                          U_star.GetLocal(global_row, column));
                }
            }
        }
    }
}
