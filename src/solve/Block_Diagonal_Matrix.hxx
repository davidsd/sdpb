//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Matrix.hxx"

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
  std::vector<Matrix> blocks;

  // The rows (or columns) of M corresponding to the top-left entry of
  // each block M_b
  std::vector<int> block_start_indices;

  // Construct a Block_Diagonal_Matrix from a vector of dimensions {s_0,
  // ..., s_{bMax-1}} for each block.
  explicit Block_Diagonal_Matrix(const std::vector<int> &block_sizes) : dim(0)
  {
    for(auto &block_size : block_sizes)
      {
        blocks.push_back(Matrix(block_size, block_size));
        block_start_indices.push_back(dim);
        dim += block_size;
      }
  }

  // M = 0
  void set_zero()
  {
    for(auto &block : blocks)
      {
        block.set_zero();
      }
  }

  // Add a constant c to each diagonal element
  void add_diagonal(const Real &c)
  {
    for(auto &block : blocks)
      {
        block.add_diagonal(c);
      }
  }

  // M = M + A
  void operator+=(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks.size(); b++)
      {
        blocks[b] += A.blocks[b];
      }
  }

  // M = M - A
  void operator-=(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks.size(); b++)
      {
        blocks[b] -= A.blocks[b];
      }
  }

  // M = c*M, where c is a constant
  void operator*=(const Real &c)
  {
    for(auto &block : blocks)
      {
        block *= c;
      }
  }

  // M = A
  void copy_from(const Block_Diagonal_Matrix &A)
  {
    for(size_t b = 0; b < blocks.size(); b++)
      {
        blocks[b].copy_from(A.blocks[b]);
      }
  }

  // Symmetrize M in place
  void symmetrize()
  {
    for(auto &block : blocks)
      {
        block.symmetrize();
      }
  }

  // The maximal absolute value of the elements of M
  Real max_abs() const
  {
    Real max = 0;
    for(auto &block : blocks)
      {
        max=std::max(block.max_abs(),max);
      }
    return max;
  }

  friend std::ostream &
  operator<<(std::ostream &os, const Block_Diagonal_Matrix &A);
};

// Tr(A B), where A and B are symmetric
Real frobeniusProductSymmetric(const Block_Diagonal_Matrix &A,
                               const Block_Diagonal_Matrix &B);

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
Real frobeniusProductOfSums(const Block_Diagonal_Matrix &X,
                            const Block_Diagonal_Matrix &dX,
                            const Block_Diagonal_Matrix &Y,
                            const Block_Diagonal_Matrix &dY);

// C := alpha*A*B + beta*C
void blockDiagonalMatrixScaleMultiplyAdd(
  Real alpha,
  Block_Diagonal_Matrix &A, // constant
  Block_Diagonal_Matrix &B, // constant
  Real beta,
  Block_Diagonal_Matrix &C); // overwritten

// C := A B
void blockDiagonalMatrixMultiply(Block_Diagonal_Matrix &A,  // constant
                                 Block_Diagonal_Matrix &B,  // constant
                                 Block_Diagonal_Matrix &C); // overwritten

// A := L^{-1} A L^{-T}
void lowerTriangularInverseCongruence(Block_Diagonal_Matrix &A,  // overwritten
                                      Block_Diagonal_Matrix &L); // constant

// Minimum eigenvalue of A, via the QR method
// Inputs:
// - A : Block_Diagonal_Matrix with blocks of size n_b x n_b (will be
//   overwritten)
// - eigenvalues : vector of Vectors of length n_b (0 <= b < bMax)
// - workspace : vector of Vectors of lenfth 3*n_b-1 (0 <= b < bMax)
//   (temporary workspace)
//
Real minEigenvalue(Block_Diagonal_Matrix &A, std::vector<Vector> &workspace,
                   std::vector<Vector> &eigenvalues);

// Compute L (lower triangular) such that A = L L^T
// Inputs:
// - A : dim x dim symmetric matrix (constant)
// - L : dim x dim lower-triangular matrix (overwritten)
//
void choleskyDecomposition(Block_Diagonal_Matrix &A, Block_Diagonal_Matrix &L);

// Compute L (lower triangular) such that A + U U^T = L L^T, where U
// is a low-rank update matrix. (See Matrix.h for details).
//
// Input:
// - A
// - stabilizeThreshold: parameter for deciding which directions to
//   stabilize
// Output:
// - L (overwritten)
// - stabilizeIndices: a list of directions which have been
//   stabilized, for each block of A (overwritten)
// - stabilizeLambdas: a list of corresponding Lambdas, for each block
//   of A (overwritten)
//
void choleskyDecompositionStabilized(
  Block_Diagonal_Matrix &A, Block_Diagonal_Matrix &L,
  std::vector<std::vector<Integer>> &stabilizeIndices,
  std::vector<std::vector<Real>> &stabilizeLambdas,
  const double stabilizeThreshold);

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
// Input:
// - ACholesky, a lower triangular Block_Diagonal_Matrix such that A =
//   ACholesky ACholesky^T (constant)
// Output:
// - X (overwritten)
void blockMatrixSolveWithCholesky(Block_Diagonal_Matrix &ACholesky, // constant
                                  Block_Diagonal_Matrix &X); // overwritten

// B := L^{-1} B, where L is lower-triangular
void blockMatrixLowerTriangularSolve(Block_Diagonal_Matrix &L, // constant
                                     Matrix &B);               // overwritten

// v := L^{-1} v, where L is lower-triangular
void blockMatrixLowerTriangularSolve(Block_Diagonal_Matrix &L, // constant
                                     Vector &v);               // overwritten

// v := L^{-T} v, where L is lower-triangular
void blockMatrixLowerTriangularTransposeSolve(
  Block_Diagonal_Matrix &L, // constant
  Vector &v);               // overwritten
