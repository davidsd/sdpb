//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include <El.hpp>

#include <list>

// A block-diagonal square matrix
//
//   M = Diagonal(M_0, M_1, ..., M_{bMax-1})
//
// where each block M_b is a square-matrix (of possibly different
// sizes).
class Block_Diagonal_Matrix
{
public:
  // The blocks M_b for 0 <= b < bMax
  std::vector<El::DistMatrix<El::BigFloat>> blocks;

  // Construct a Block_Diagonal_Matrix from a vector of dimensions {s_0,
  // ..., s_{bMax-1}} for each block.
  explicit Block_Diagonal_Matrix(const std::vector<size_t> &block_sizes,
                                 const std::vector<size_t> &block_indices,
                                 const size_t &num_schur_blocks,
                                 const El::Grid &grid)
  {
    bool scale_index(num_schur_blocks != block_sizes.size());
    blocks.reserve(block_indices.size() * (scale_index ? 2 : 1));
    for(auto &block_index : block_indices)
      {
        if(scale_index)
          {
            add_block(block_sizes.at(block_index * 2), grid);
            add_block(block_sizes.at(block_index * 2 + 1), grid);
          }
        else
          {
            add_block(block_sizes.at(block_index), grid);
          }
      }
  }

  void add_block(const size_t &block_size, const El::Grid &grid)
  {
    blocks.emplace_back(block_size, block_size, grid);
  }

  void set_zero()
  {
    for(auto &block : blocks)
      {
        El::Zero(block);
      }
  }

  // Add a constant c to each diagonal element
  void add_diagonal(const El::BigFloat &c)
  {
    for(auto &block : blocks)
      {
        ShiftDiagonal(block, c);
      }
  }

  void operator+=(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks.size(); b++)
      {
        blocks[b] += A.blocks[b];
      }
  }

  void operator-=(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks.size(); b++)
      {
        blocks[b] -= A.blocks[b];
      }
  }

  void operator*=(const El::BigFloat &c)
  {
    for(auto &block : blocks)
      {
        block *= c;
      }
  }

  void symmetrize()
  {
    for(auto &block : blocks)
      {
        // FIXME: This feels expensive

        // We can not use El::MakeSymmetric() because that just copies
        // the lower part to the upper part.  We need to average the
        // upper and lower parts.
        block *= 0.5;
        El::DistMatrix<El::BigFloat> transpose(block.Grid());
        El::Transpose(block, transpose, false);
        block += transpose;
      }
  }

  // The maximal absolute value of the elements of M
  El::BigFloat max_abs() const
  {
    El::BigFloat max = 0;
    for(auto &block : blocks)
      {
        max = std::max(El::MaxAbs(block), max);
      }
    return max;
  }

  El::BigFloat max_abs_mpi() const
  {
	  El::BigFloat local_max = 0;
	  for (auto &block : blocks)
	  {
		  local_max = std::max(El::MaxAbs(block), local_max);
	  }
	  return El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);
  }

  friend std::ostream &
  operator<<(std::ostream &os, const Block_Diagonal_Matrix &A);
};

