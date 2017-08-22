//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_BLOCKDIAGONALMATRIX_H_
#define SDPB_BLOCKDIAGONALMATRIX_H_

#include <vector>
#include <iostream>
#include <ostream>
#include "omp.h"
#include "types.h"
#include "Matrix.h"

using std::vector;
using std::ostream;

// A block-diagonal square matrix
//
//   M = Diagonal(M_0, M_1, ..., M_{bMax-1})
//
// where each block M_b is a square-matrix (of possibly different
// sizes).
class BlockDiagonalMatrix {
public:
  // The total number of rows (or columns) in M
  int dim;

  // The blocks M_b for 0 <= b < bMax
  vector<Matrix> blocks;

  // The rows (or columns) of M corresponding to the top-left entry of
  // each block M_b
  vector<int> blockStartIndices;

  // Construct a BlockDiagonalMatrix from a vector of dimensions {s_0,
  // ..., s_{bMax-1}} for each block.
  explicit BlockDiagonalMatrix(const vector<int> &blockSizes):
    dim(0) {
    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(Matrix(blockSizes[i], blockSizes[i]));
      blockStartIndices.push_back(dim);
      dim += blockSizes[i];
    }
  }

  // M = 0
  void setZero() {
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].setZero();
  }

  // Add a constant c to each diagonal element
  void addDiagonal(const Real &c) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].addDiagonal(c);
  }

  // M = M + A
  void operator+=(const BlockDiagonalMatrix &A) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] += A.blocks[b];
  }

  // M = M - A
  void operator-=(const BlockDiagonalMatrix &A) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] -= A.blocks[b];
  }

  // M = c*M, where c is a constant
  void operator*=(const Real &c) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] *= c;
  }

  // M = A
  void copyFrom(const BlockDiagonalMatrix &A) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].copyFrom(A.blocks[b]);
  }

  // Symmetrize M in place
  void symmetrize() {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].symmetrize();
  }

  // The maximal absolute value of the elements of M
  Real maxAbs() const {
    Real max = 0;
    for (vector<Matrix>::const_iterator b = blocks.begin();
         b != blocks.end();
         b++) {
      const Real tmp = b->maxAbs();
      if (tmp > max)
        max = tmp;
    }
    return max;
  }

  friend ostream& operator<<(ostream& os, const BlockDiagonalMatrix& A);
};

// Tr(A B), where A and B are symmetric
Real frobeniusProductSymmetric(const BlockDiagonalMatrix &A,
                               const BlockDiagonalMatrix &B);

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
Real frobeniusProductOfSums(const BlockDiagonalMatrix &X,
                            const BlockDiagonalMatrix &dX,
                            const BlockDiagonalMatrix &Y,
                            const BlockDiagonalMatrix &dY);

// C := alpha*A*B + beta*C
void blockDiagonalMatrixScaleMultiplyAdd(Real alpha,
                                         BlockDiagonalMatrix &A,  // constant
                                         BlockDiagonalMatrix &B,  // constant 
                                         Real beta,
                                         BlockDiagonalMatrix &C); // overwritten

// C := A B
void blockDiagonalMatrixMultiply(BlockDiagonalMatrix &A,  // constant
                                 BlockDiagonalMatrix &B,  // constant
                                 BlockDiagonalMatrix &C); // overwritten

// A := L^{-1} A L^{-T}
void lowerTriangularInverseCongruence(BlockDiagonalMatrix &A,  // overwritten
                                      BlockDiagonalMatrix &L); // constant

// Minimum eigenvalue of A, via the QR method
// Inputs:
// - A : BlockDiagonalMatrix with blocks of size n_b x n_b (will be
//   overwritten)
// - eigenvalues : vector of Vectors of length n_b (0 <= b < bMax)
// - workspace : vector of Vectors of lenfth 3*n_b-1 (0 <= b < bMax)
//   (temporary workspace)
//
Real minEigenvalue(BlockDiagonalMatrix &A,
                   vector<Vector> &workspace,
                   vector<Vector> &eigenvalues);

// Compute L (lower triangular) such that A = L L^T
// Inputs:
// - A : dim x dim symmetric matrix (constant)
// - L : dim x dim lower-triangular matrix (overwritten)
//
void choleskyDecomposition(BlockDiagonalMatrix &A,
                           BlockDiagonalMatrix &L);

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
void choleskyDecompositionStabilized(BlockDiagonalMatrix &A,
                                     BlockDiagonalMatrix &L,
                                     vector<vector<Integer> > &stabilizeIndices,
                                     vector<vector<Real> > &stabilizeLambdas,
                                     const double stabilizeThreshold);

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
// Input:
// - ACholesky, a lower triangular BlockDiagonalMatrix such that A =
//   ACholesky ACholesky^T (constant)
// Output:
// - X (overwritten)
void blockMatrixSolveWithCholesky(BlockDiagonalMatrix &ACholesky, // constant
                                  BlockDiagonalMatrix &X);        // overwritten

// B := L^{-1} B, where L is lower-triangular
void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, // constant
                                     Matrix &B);             // overwritten

// v := L^{-1} v, where L is lower-triangular
void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, // constant
                                     Vector &v);             // overwritten

// v := L^{-T} v, where L is lower-triangular
void blockMatrixLowerTriangularTransposeSolve(BlockDiagonalMatrix &L, // constant
                                              Vector &v);             // overwritten

#endif  // SDPB_BLOCKDIAGONALMATRIX_H_
