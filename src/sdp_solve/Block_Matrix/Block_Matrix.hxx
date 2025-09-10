#pragma once

// A general matrix with the same vertical distribution as
// Block_Diagonal_Matrix
//
//    ( M_0 )
// M =( M_1 )
//    ( M_2 )
//    ( ... )
//
// The blocks are not, in general, square.  This allows us to compute
// solutions for each block independently.

#include "Abstract_Block_Matrix.hxx"

struct Block_Matrix : Abstract_Block_Matrix<Block_Matrix>
{
  virtual ~Block_Matrix() = default;
  size_t width;
  Block_Matrix(const std::vector<size_t> &block_heights, const size_t &width,
               const std::vector<size_t> &block_indices, const El::Grid &grid)
      : width(width)
  {
    blocks.reserve(block_indices.size());
    for(auto &block_index : block_indices)
      {
        Block_Matrix::add_block(block_heights[block_index], grid);
      }
  }

  void add_block(size_t height, const El::Grid &grid) override
  {
    blocks.emplace_back(height, width, grid);
  }

  Block_Matrix() = default;
};
