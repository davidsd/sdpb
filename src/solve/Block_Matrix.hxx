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
               const std::list<El::Grid> &block_grid_mapping)
  {
    bool same_size(block_grid_mapping.size() == block_heights.size()),
      skip_increment(true);
    auto grid(block_grid_mapping.begin());
    for(auto &block_height : block_heights)
      {
        blocks.emplace_back(block_height, width, *grid);
        if(same_size)
          {
            ++grid;
          }
        else if(!same_size)
          {
            if(!skip_increment)
              {
                ++grid;
              }
            skip_increment = !skip_increment;
          }
      }
  }
  Block_Matrix() = default;

  size_t width() const
  {
    size_t result(0);
    if(!blocks.empty())
      {
        result = blocks[0].Width();
      }
    return result;
  }
};
