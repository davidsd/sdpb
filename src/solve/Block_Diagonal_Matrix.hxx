//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Matrix.hxx"
#include <El.hpp>

// A block-diagonal square matrix
//
//   M = Diagonal(M_0, M_1, ..., M_{bMax-1})
//
// where each block M_b is a square-matrix (of possibly different
// sizes).
class Block_Diagonal_Matrix
{
public:
  // The total number of rows (or columns) in M
  int dim;

  // The blocks M_b for 0 <= b < bMax
  std::vector<El::DistMatrix<El::BigFloat>> blocks_elemental;

  // The rows (or columns) of M corresponding to the top-left entry of
  // each block M_b
  std::vector<size_t> block_start_indices;

  // Construct a Block_Diagonal_Matrix from a vector of dimensions {s_0,
  // ..., s_{bMax-1}} for each block.
  explicit Block_Diagonal_Matrix(const std::vector<size_t> &block_sizes)
      : dim(0)
  {
    for(auto &block_size : block_sizes)
      {
        blocks_elemental.emplace_back(block_size, block_size);

        block_start_indices.push_back(dim);
        dim += block_size;
      }
  }

  // M = 0
  void set_zero()
  {
    for(auto &block : blocks_elemental)
      {
        Zero(block);
      }
  }

  // Add a constant c to each diagonal element
  void add_diagonal(const El::BigFloat &c)
  {
    for(auto &block : blocks_elemental)
      {
        ShiftDiagonal(block, c);
      }
  }

  // M = M + A
  void operator+=(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks_elemental.size(); b++)
      {
        blocks_elemental[b] += A.blocks_elemental[b];
      }
  }

  // M = M - A
  void operator-=(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks_elemental.size(); b++)
      {
        blocks_elemental[b] -= A.blocks_elemental[b];
      }
  }

  // M = c*M, where c is a constant
  void operator*=(const El::BigFloat &c)
  {
    for(auto &block : blocks_elemental)
      {
        block *= c;
      }
  }

  // M = A
  void copy_from(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks_elemental.size(); b++)
      {
        El::Copy(A.blocks_elemental[b], blocks_elemental[b]);
      }
  }

  // Symmetrize M in place
  void symmetrize()
  {
    for(auto &block : blocks_elemental)
      {
        // FIXME: This feels expensive

        // We can not use El::MakeSymmetric() because that just copies
        // the lower part to the upper part.  We need to average the
        // upper and lower parts.
        block *= 0.5;
        El::DistMatrix<El::BigFloat> transpose;
        El::Transpose(block, transpose, false);
        block += transpose;
      }
  }

  // The maximal absolute value of the elements of M
  El::BigFloat max_abs() const
  {
    El::BigFloat max = 0;
    for(auto &block : blocks_elemental)
      {
        max = std::max(El::MaxAbs(block), max);
      }
    return max;
  }

  friend std::ostream &
  operator<<(std::ostream &os, const Block_Diagonal_Matrix &A);
};

// Tr(A B), where A and B are symmetric
El::BigFloat
frobenius_product_symmetric_elemental(const Block_Diagonal_Matrix &A,
                                      const Block_Diagonal_Matrix &B);

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
El::BigFloat
frobenius_product_of_sums_elemental(const Block_Diagonal_Matrix &X,
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
void block_diagonal_matrix_multiply(Block_Diagonal_Matrix &A,  // constant
                                    Block_Diagonal_Matrix &B,  // constant
                                    Block_Diagonal_Matrix &C); // overwritten

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(
  Block_Diagonal_Matrix &A,  // overwritten
  Block_Diagonal_Matrix &L); // constant

// Minimum eigenvalue of A, via the QR method
// Inputs:
// - A : Block_Diagonal_Matrix with blocks of size n_b x n_b (will be
//   overwritten)
// - eigenvalues : vector of Vectors of length n_b (0 <= b < bMax)
// - workspace : vector of Vectors of lenfth 3*n_b-1 (0 <= b < bMax)
//   (temporary workspace)
//
El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A);

// Compute L (lower triangular) such that A = L L^T
// Inputs:
// - A : dim x dim symmetric matrix (constant)
// - L : dim x dim lower-triangular matrix (overwritten)
//
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
// Input:
// - ACholesky, a lower triangular Block_Diagonal_Matrix such that A =
//   ACholesky ACholesky^T (constant)
// Output:
// - X (overwritten)
void block_matrix_solve_with_cholesky(
  const Block_Diagonal_Matrix &ACholesky, // constant
  Block_Diagonal_Matrix &X);              // overwritten
