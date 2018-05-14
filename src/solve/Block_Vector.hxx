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

#include <El.hpp>

#include <vector>

struct Block_Vector
{
  std::vector<El::DistMatrix<El::BigFloat>> blocks;

  Block_Vector(const std::vector<size_t> &block_heights)
  {
    for(auto &block_height : block_heights)
      {
        blocks.emplace_back(block_height, 1);
      }
  }
  Block_Vector() = default;
};
