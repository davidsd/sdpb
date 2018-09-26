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
        Zero(block);
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

  friend std::ostream &
  operator<<(std::ostream &os, const Block_Diagonal_Matrix &A);
};

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
El::BigFloat frobenius_product_of_sums(const Block_Diagonal_Matrix &X,
                                       const Block_Diagonal_Matrix &dX,
                                       const Block_Diagonal_Matrix &Y,
                                       const Block_Diagonal_Matrix &dY);

// C := alpha*A*B + beta*C
void block_diagonal_matrix_scale_multiply_add(const El::BigFloat &alpha,
                                              const Block_Diagonal_Matrix &A,
                                              const Block_Diagonal_Matrix &B,
                                              const El::BigFloat &beta,
                                              Block_Diagonal_Matrix &C);

// C := A B
void block_diagonal_matrix_multiply(const Block_Diagonal_Matrix &A,
                                    const Block_Diagonal_Matrix &B,
                                    Block_Diagonal_Matrix &C);

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(const Block_Diagonal_Matrix &L,
                                         Block_Diagonal_Matrix &A);

El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A);

void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void block_matrix_solve_with_cholesky(
  const Block_Diagonal_Matrix &ACholesky,
  Block_Diagonal_Matrix &X);
