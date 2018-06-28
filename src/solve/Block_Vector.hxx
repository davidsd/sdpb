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

#include <list>
#include <vector>

struct Block_Vector
{
  std::vector<El::DistMatrix<El::BigFloat>> blocks;

  Block_Vector(const std::vector<size_t> &block_heights,
               const std::list<El::Grid> &block_grid_mapping)
  {
    bool same_size(block_grid_mapping.size() == block_heights.size()),
      skip_increment(true);
    auto grid(block_grid_mapping.begin());
    for(auto &block_height : block_heights)
      {
        blocks.emplace_back(block_height, 1, *grid);
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
  Block_Vector() = default;
};
