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

#include <vector>

struct Block_Matrix
{
  std::vector<El::DistMatrix<El::BigFloat>> blocks_elemental;

  Block_Matrix(const std::vector<int> &block_heights, const size_t &width)
  {
    for(auto &block_height : block_heights)
      {
        blocks_elemental.emplace_back(block_height, width);
      }
  }
};
