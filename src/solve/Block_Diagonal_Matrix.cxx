//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "Block_Diagonal_Matrix.hxx"
#include <algorithm>
#include <vector>

std::ostream &operator<<(std::ostream &os, const Block_Diagonal_Matrix &A)
{
  os << "{";
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    {
      os << A.blocks[b];
      if(b < A.blocks.size() - 1)
        os << ", ";
    }
  os << "}";
  return os;
}

// Tr(A B), where A and B are symmetric
Real frobeniusProductSymmetric(const Block_Diagonal_Matrix &A,
                               const Block_Diagonal_Matrix &B)
{
  Real result = 0;
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    {
      Real f = frobeniusProductSymmetric(A.blocks[b], B.blocks[b]);
      // this pragma means that other threads should be stopped while
      // this operation is performed (this prevents result from being
      // modified by multiple threads simultaneously)
      {
        result += f;
      }
    }
  return result;
}

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
Real frobeniusProductOfSums(const Block_Diagonal_Matrix &X,
                            const Block_Diagonal_Matrix &dX,
                            const Block_Diagonal_Matrix &Y,
                            const Block_Diagonal_Matrix &dY)
{
  Real result = 0;
  for(unsigned int b = 0; b < X.blocks.size(); b++)
    {
      Real f = frobeniusProductOfSums(X.blocks[b], dX.blocks[b], Y.blocks[b],
                                      dY.blocks[b]);
      {
        result += f;
      }
    }
  return result;
}

// C := alpha*A*B + beta*C
void blockDiagonalMatrixScaleMultiplyAdd(Real alpha, Block_Diagonal_Matrix &A,
                                         Block_Diagonal_Matrix &B, Real beta,
                                         Block_Diagonal_Matrix &C)
{
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    matrixScaleMultiplyAdd(alpha, A.blocks[b], B.blocks[b], beta, C.blocks[b]);
}

// C := A*B
void blockDiagonalMatrixMultiply(Block_Diagonal_Matrix &A,
                                 Block_Diagonal_Matrix &B,
                                 Block_Diagonal_Matrix &C)
{
  blockDiagonalMatrixScaleMultiplyAdd(1, A, B, 0, C);
}

// A := L^{-1} A L^{-T}
void lowerTriangularInverseCongruence(Block_Diagonal_Matrix &A,
                                      Block_Diagonal_Matrix &L)
{
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    lowerTriangularInverseCongruence(A.blocks[b], L.blocks[b]);
}

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric Block_Diagonal_Matrix
// eigenvalues : vector<Vector> of length A.blocks.size()
// workspace   : vector<Vector> of length A.blocks.size()
//
Real minEigenvalue(Block_Diagonal_Matrix &A, vector<Vector> &workspace,
                   vector<Vector> &eigenvalues)
{
  assert(A.blocks.size() == eigenvalues.size());
  assert(A.blocks.size() == workspace.size());

  // TODO(davidsd): get rid of this hack
  Real lambdaMin = 1e100; // we really want lambdaMin = infinity here
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    {
      Real minBlockLambda
        = minEigenvalue(A.blocks[b], workspace[b], eigenvalues[b]);
      // ensure only one thread modifies lambdaMin at a time
      {
        lambdaMin = min(lambdaMin, minBlockLambda);
      }
    }

  return lambdaMin;
}

// Compute L (lower triangular) such that A = L L^T
void choleskyDecomposition(Block_Diagonal_Matrix &A, Block_Diagonal_Matrix &L)
{
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    choleskyDecomposition(A.blocks[b], L.blocks[b]);
}

// Compute L (lower triangular) such that A + U U^T = L L^T, where U
// is a low-rank update matrix.
void choleskyDecompositionStabilized(Block_Diagonal_Matrix &A,
                                     Block_Diagonal_Matrix &L,
                                     vector<vector<Integer>> &stabilizeIndices,
                                     vector<vector<Real>> &stabilizeLambdas,
                                     const double stabilizeThreshold)
{
  for(unsigned int b = 0; b < A.blocks.size(); b++)
    choleskyDecompositionStabilized(A.blocks[b], L.blocks[b],
                                    stabilizeIndices[b], stabilizeLambdas[b],
                                    stabilizeThreshold);
}

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void blockMatrixSolveWithCholesky(Block_Diagonal_Matrix &ACholesky,
                                  Block_Diagonal_Matrix &X)
{
  for(unsigned int b = 0; b < X.blocks.size(); b++)
    matrixSolveWithCholesky(ACholesky.blocks[b], X.blocks[b]);
}

// B := L^{-1} B, where L is lower-triangular
void blockMatrixLowerTriangularSolve(Block_Diagonal_Matrix &L, Matrix &B)
{
  for(unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularSolve(L.blocks[b], &B.elt(L.blockStartIndices[b], 0),
                         B.cols, B.rows);
}

// v := L^{-1} v, where L is lower-triangular
void blockMatrixLowerTriangularSolve(Block_Diagonal_Matrix &L, Vector &v)
{
  for(unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularSolve(L.blocks[b], &v[L.blockStartIndices[b]], 1, v.size());
}

// v := L^{-T} v, where L is lower-triangular
void blockMatrixLowerTriangularTransposeSolve(Block_Diagonal_Matrix &L,
                                              Vector &v)
{
  for(unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularTransposeSolve(L.blocks[b], &v[L.blockStartIndices[b]], 1,
                                  v.size());
}
