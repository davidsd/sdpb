#include "BlockDiagonalMatrix.h"

ostream& operator<<(ostream& os, const BlockDiagonalMatrix& A) {
  os << "{";
  for (unsigned int b = 0; b < A.blocks.size(); b++) {
    os << A.blocks[b];
    if (b < A.blocks.size() - 1)
      os << ", ";
  }
  os << "}";
  return os;
}

Real frobeniusProductSymmetric(const BlockDiagonalMatrix &A,
                               const BlockDiagonalMatrix &B) {
  Real result = 0;
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++) {
    Real f = frobeniusProductSymmetric(A.blocks[b], B.blocks[b]);
    #pragma omp critical
    {
      result += f;
    }
  }
  return result;
}
  
// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
Real frobeniusProductOfSums(const BlockDiagonalMatrix &X,
                            const BlockDiagonalMatrix &dX,
                            const BlockDiagonalMatrix &Y,
                            const BlockDiagonalMatrix &dY) {
  Real result = 0;
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < X.blocks.size(); b++) {
    Real f = frobeniusProductOfSums(X.blocks[b], dX.blocks[b], Y.blocks[b], dY.blocks[b]);
    #pragma omp critical
    {
      result += f;
    }
  }
  return result;
}

void blockDiagonalMatrixScaleMultiplyAdd(Real alpha,
                                         BlockDiagonalMatrix &A,
                                         BlockDiagonalMatrix &B,
                                         Real beta,
                                         BlockDiagonalMatrix &C) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    matrixScaleMultiplyAdd(alpha, A.blocks[b], B.blocks[b], beta, C.blocks[b]);
}

void blockDiagonalMatrixMultiply(BlockDiagonalMatrix &A,
                                 BlockDiagonalMatrix &B,
                                 BlockDiagonalMatrix &C) {
  blockDiagonalMatrixScaleMultiplyAdd(1, A, B, 0, C);
}

void lowerTriangularInverseCongruence(BlockDiagonalMatrix &A, BlockDiagonalMatrix &L) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    lowerTriangularInverseCongruence(A.blocks[b], L.blocks[b]);
}

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric BlockDiagonalMatrix
// eigenvalues : vector<Vector> of length A.blocks.size()
// workspace   : vector<Vector> of length A.blocks.size()
//
Real minEigenvalue(BlockDiagonalMatrix &A, vector<Vector> &workspace, vector<Vector> &eigenvalues) {
  assert(A.blocks.size() == eigenvalues.size());
  assert(A.blocks.size() == workspace.size());

  // TODO: get rid of this hack
  Real lambdaMin = 1e100;
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++) {
    Real minBlockLambda = minEigenvalue(A.blocks[b], workspace[b], eigenvalues[b]);
    #pragma omp critical
    {
      lambdaMin = min(lambdaMin, minBlockLambda);
    }
  }

  return lambdaMin;
}

void choleskyDecomposition(BlockDiagonalMatrix &A,
                           BlockDiagonalMatrix &L) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    choleskyDecomposition(A.blocks[b], L.blocks[b]);
}

void choleskyDecompositionStabilized(BlockDiagonalMatrix &A,
                                     BlockDiagonalMatrix &L,
                                     vector<vector<Integer> > &schurStabilizeIndices,
                                     vector<vector<Real> > &schurStabilizeLambdas) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    choleskyDecompositionStabilized(A.blocks[b], L.blocks[b], schurStabilizeIndices[b], schurStabilizeLambdas[b]);
}

void blockMatrixSolveWithCholesky(BlockDiagonalMatrix &ACholesky,
                                  BlockDiagonalMatrix &X) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < X.blocks.size(); b++)
    matrixSolveWithCholesky(ACholesky.blocks[b], X.blocks[b]);
}

void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, Matrix &B) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularSolve(L.blocks[b], &B.elt(L.blockStartIndices[b], 0), B.cols, B.rows);
}

void blockMatrixLowerTriangularTransposeSolve(BlockDiagonalMatrix &L, Matrix &B) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularTransposeSolve(L.blocks[b], &B.elt(L.blockStartIndices[b], 0), B.cols, B.rows);
}

void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, Vector &v) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularSolve(L.blocks[b], &v[L.blockStartIndices[b]], 1, v.size());
}

void blockMatrixLowerTriangularTransposeSolve(BlockDiagonalMatrix &L, Vector &v) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularTransposeSolve(L.blocks[b], &v[L.blockStartIndices[b]], 1, v.size());
}
