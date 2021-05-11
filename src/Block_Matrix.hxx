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

#include <El.hpp>

#include <list>
#include <vector>

struct Block_Matrix
{
  std::vector<El::DistMatrix<El::BigFloat>> blocks;

  Block_Matrix(const std::vector<size_t> &block_heights, const size_t &width,
               const std::vector<size_t> &block_indices,
               const size_t &num_schur_blocks, const El::Grid &grid)
  {
    bool scale_index(num_schur_blocks != block_heights.size());
    blocks.reserve(block_indices.size() * (scale_index ? 2 : 1));
    for(auto &block_index : block_indices)
      {
        if(scale_index)
          {
            blocks.emplace_back(block_heights.at(2 * block_index), width,
                                grid);
            blocks.emplace_back(block_heights.at(2 * block_index + 1), width,
                                grid);
          }
        else
          {
            blocks.emplace_back(block_heights[block_index], width, grid);
          }
      }
  }
  Block_Matrix() = default;
};
