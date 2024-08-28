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

template <class TDistMatrix>
struct Block_Vector_Impl
{
  std::vector<TDistMatrix> blocks;

  Block_Vector_Impl(const std::vector<size_t> &block_heights,
                    const std::vector<size_t> &block_indices,
                    const size_t &num_schur_blocks, const El::Grid &grid)
  {
    bool scale_index(num_schur_blocks != block_heights.size());
    blocks.reserve(block_indices.size() * (scale_index ? 2 : 1));
    for(auto &block_index : block_indices)
      {
        if(scale_index)
          {
            blocks.emplace_back(block_heights.at(2 * block_index), 1, grid);
            blocks.emplace_back(block_heights.at(2 * block_index + 1), 1,
                                grid);
          }
        else
          {
            blocks.emplace_back(block_heights.at(block_index), 1, grid);
          }
      }
  }
  Block_Vector_Impl() = default;
};

using Block_Vector = Block_Vector_Impl<El::DistMatrix<El::BigFloat>>;

// Each element is copied over all ranks in group
using Block_Vector_Star
  = Block_Vector_Impl<El::DistMatrix<El::BigFloat, El::STAR, El::STAR>>;
