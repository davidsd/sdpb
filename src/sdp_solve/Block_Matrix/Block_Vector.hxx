#pragma once

// A general vector with the same vertical distribution as
// Block_Diagonal_Matrix
//
//    ( M_0 )
// M =( M_1 )
//    ( M_2 )
//    ( ... )
//
// This is equivalent to Block_Matrix with width=1.  We use a separate
// type to enhance type safety.

#include "Abstract_Block_Matrix.hxx"

struct Block_Vector : Abstract_Block_Matrix<Block_Vector>
{
  Block_Vector(const std::vector<size_t> &block_heights,
               const std::vector<size_t> &block_indices, const El::Grid &grid)
  {
    blocks.reserve(block_indices.size());
    for(auto &block_index : block_indices)
      {
        Block_Vector::add_block(block_heights.at(block_index), grid);
      }
  }
  void add_block(size_t height, const El::Grid &grid) override
  {
    blocks.emplace_back(height, 1, grid);
  }
  Block_Vector() = default;
};
