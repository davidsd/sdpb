//=======================================================================
// Copyright 2014 David Simmons-Duffin.
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

class BlockDiagonalMatrix {
public:
  int dim;
  vector<Matrix> blocks;
  vector<int> blockStartIndices;

  BlockDiagonalMatrix(const vector<int> &blockSizes):
    dim(0) {
    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(Matrix(blockSizes[i], blockSizes[i]));
      blockStartIndices.push_back(dim);
      dim += blockSizes[i];
    }
  }

  void setZero() {
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].setZero();
  }

  void addDiagonal(const Real &c) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].addDiagonal(c);
  }

  void operator+=(const BlockDiagonalMatrix &A) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] += A.blocks[b];
  }

  void operator-=(const BlockDiagonalMatrix &A) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] -= A.blocks[b];
  }

  void operator*=(const Real &c) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] *= c;
  }

  void copyFrom(const BlockDiagonalMatrix &A) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].copyFrom(A.blocks[b]);
  }

  void symmetrize() {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].symmetrize();
  }

  Real maxAbs() const {
    Real max = 0;
    for (vector<Matrix>::const_iterator b = blocks.begin(); b != blocks.end(); b++) {
      const Real tmp = b->maxAbs();
      if (tmp > max)
        max = tmp;
    }
    return max;
  }

  friend ostream& operator<<(ostream& os, const BlockDiagonalMatrix& A);

};

Real frobeniusProductSymmetric(const BlockDiagonalMatrix &A,
                               const BlockDiagonalMatrix &B);
  
// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
Real frobeniusProductOfSums(const BlockDiagonalMatrix &X,
                            const BlockDiagonalMatrix &dX,
                            const BlockDiagonalMatrix &Y,
                            const BlockDiagonalMatrix &dY);

void blockDiagonalMatrixScaleMultiplyAdd(Real alpha,
                                         BlockDiagonalMatrix &A,
                                         BlockDiagonalMatrix &B,
                                         Real beta,
                                         BlockDiagonalMatrix &C);
void blockDiagonalMatrixMultiply(BlockDiagonalMatrix &A, BlockDiagonalMatrix &B, BlockDiagonalMatrix &C);
void lowerTriangularInverseCongruence(BlockDiagonalMatrix &A, BlockDiagonalMatrix &L);
Real minEigenvalue(BlockDiagonalMatrix &A, vector<Vector> &workspace, vector<Vector> &eigenvalues);
void choleskyDecomposition(BlockDiagonalMatrix &A, BlockDiagonalMatrix &L);
void choleskyDecompositionStabilized(BlockDiagonalMatrix &A,
                                     BlockDiagonalMatrix &L,
                                     vector<vector<Integer> > &schurStabilizeIndices,
                                     vector<vector<Real> > &schurStabilizeLambdas,
                                     const double stabilizeThreshold);
void blockMatrixSolveWithCholesky(BlockDiagonalMatrix &ACholesky, BlockDiagonalMatrix &X);
void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, Matrix &B);
void blockMatrixLowerTriangularTransposeSolve(BlockDiagonalMatrix &L, Matrix &B);
void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, Vector &v);
void blockMatrixLowerTriangularTransposeSolve(BlockDiagonalMatrix &L, Vector &v);

#endif  // SDPB_MATRIX_H_
