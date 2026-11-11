//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Abstract_Block_Diagonal_Matrix.hxx"
#include "sdpb_util/assert.hxx"

// A block-diagonal square matrix
//
//   M = Diagonal(M_0, M_1, ..., M_{bMax-1})
//
// where each block M_b is a square-matrix (of possibly different
// sizes).
struct Block_Diagonal_Matrix : Abstract_Block_Diagonal_Matrix<Block_Diagonal_Matrix>
{
  // Construct a Block_Diagonal_Matrix from a vector of dimensions {s_0,
  // ..., s_{bMax-1}} for each block.
  Block_Diagonal_Matrix(const std::vector<size_t> &block_sizes,
                        const std::vector<size_t> &block_indices,
                        const El::Grid &grid);
};

// A block-diagonal distributed square matrix, where blocks are grouped in pairs.
// M = Diagonal(M_00, M_01, M_10, M_11, ..., M_{bMax-1}0, M_{bMax-1}1)
// where each block M_bi is a square-matrix (of possibly different sizes).
// M_b0 and M_b1 have the same size.
struct Paired_Block_Diagonal_Matrix : Abstract_Block_Diagonal_Matrix<Paired_Block_Diagonal_Matrix>
{
  // Construct a Paired_Block_Diagonal_Matrix from a vector of dimensions {s_0,
  // ..., s_{bMax-1}} for each block.
  // Each block_index from block_indices gets two blocks in the matrix.
  Paired_Block_Diagonal_Matrix(const std::vector<size_t> &block_sizes,
                               const std::vector<size_t> &block_indices,
                               const El::Grid &grid);

  El::DistMatrix<El::BigFloat> &get_block(size_t block_index, size_t parity);

  const El::DistMatrix<El::BigFloat> &
  get_block(size_t block_index, size_t parity) const;
};

// Compute L (lower triangular) such that A = L L^T
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L,
                            const std::vector<size_t> &block_indices,
                            const std::string &name);

// Compute L (lower triangular) such that A = L L^T
void cholesky_decomposition(const Paired_Block_Diagonal_Matrix &A,
                            Paired_Block_Diagonal_Matrix &L,
                            const std::vector<size_t> &block_indices,
                            const std::string &name);
