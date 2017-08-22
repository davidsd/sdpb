//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_MATRIX_H_
#define SDPB_MATRIX_H_

#include <assert.h>
#include <iostream>
#include <ostream>
#include <vector>
#include "omp.h"
#include "types.h"
#include "Vector.h"

using std::ostream;
using std::vector;

// A matrix M with Real entries
class Matrix {
 public:
  int rows;
  int cols;
  // Elements of M in row-major order:
  //
  //   elements = { M_{0,0}, ..., M_{cols-1,0},
  //                M_{0,1}, ..., M_{cols-1,1},
  //                ...
  //                M_{0,rows-1}, ..., M_{cols-1,rows-1} }
  //
  Vector elements;

  Matrix(int rows = 0, int cols = 0):
    rows(rows),
    cols(cols),
    elements(Vector(rows*cols, 0)) {}

  inline const Real& elt(const int r, const int c) const {
    return elements[r + c*rows];
  }

  inline Real& elt(const int r, const int c) {
    return elements[r + c*rows];
  }

  // M := 0
  void setZero() {
    std::fill(elements.begin(), elements.end(), 0);
  }

  // M += c*I, where I is the identity and c is a constant
  void addDiagonal(const Real &c) {
    // ensure M is square
    assert(rows == cols);

    for (int i = 0; i < rows; i++)
      elt(i, i) += c;
  }

  void resize(int r, int c) {
    elements.resize(r*c);
    rows = r;
    cols = c;
  }

  void symmetrize() {
    assert(rows == cols);

    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < r; c++) {
        Real tmp = (elt(r, c) + elt(c, r))/2;
        elt(r, c) = tmp;
        elt(c, r) = tmp;
      }
    }
  }

  // M := A
  void copyFrom(const Matrix &A) {
    assert(rows == A.rows);
    assert(cols == A.cols);

    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] = A.elements[i];
  }

  // M := M + A
  void operator+=(const Matrix &A) {
    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] += A.elements[i];
  }

  // M := M - A
  void operator-=(const Matrix &A) {
    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] -= A.elements[i];
  }

  // M := c*M, where c is a constant
  void operator*=(const Real &c) {
    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] *= c;
  }

  // The maximum absolute value of the elemnts of M
  Real maxAbs() const {
    return maxAbsVector(elements);
  }

  friend ostream& operator<<(ostream& os, const Matrix& a);
};

// C := alpha*A*B + beta*C
void matrixScaleMultiplyAdd(Real alpha, Matrix &A, Matrix &B,
                            Real beta, Matrix &C);

// C := A*B
void matrixMultiply(Matrix &A, Matrix &B, Matrix &C);

// Set block starting at (bRow, bCol) of B to A^T A
void matrixSquareIntoBlock(Matrix &A, Matrix &B, int bRow, int bCol);

// A := L^{-1} A L^{-T}
void lowerTriangularInverseCongruence(Matrix &A, Matrix &L);

// y := alpha A x + beta y
void vectorScaleMatrixMultiplyAdd(Real alpha, Matrix &A, Vector &x,
                                  Real beta, Vector &y);

// y := alpha A^T x + beta y
void vectorScaleMatrixMultiplyTransposeAdd(Real alpha, Matrix &A, Vector &x,
                                           Real beta, Vector &y);

// Frobenius product Tr(A^T B) where A and B are symmetric matrices
Real frobeniusProductSymmetric(const Matrix &A, const Matrix &B);

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric Matrices and
// '.' is the Frobenius product.
Real frobeniusProductOfSums(const Matrix &X, const Matrix &dX,
                            const Matrix &Y, const Matrix &dY);

// Eigenvalues of A, via the QR method
// Inputs:
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
void matrixEigenvalues(Matrix &A, Vector &workspace, Vector &eigenvalues);

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : n x n Matrix (overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
Real minEigenvalue(Matrix &A, Vector &workspace, Vector &eigenvalues);

// Compute an in-place LU decomposition of A, with pivots, suitable
// for use with 'solveWithLUDecomposition'
void LUDecomposition(Matrix &A, vector<Integer> &pivots);

// b := A^{-1} b, where LU and pivots encode the LU decomposition of A
void solveWithLUDecomposition(Matrix &LU, vector<Integer> &pivots, Vector &b);

// L (lower triangular) such that A = L L^T
// Input:
// - A : dim x dim symmetric matrix
// Output:
// - L : dim x dim lower-triangular matrix (overwritten)
void choleskyDecomposition(Matrix &A, Matrix &L);

// Compute L (lower triangular) such that A + U U^T = L L^T.  Here,
// the 'update' matrix U has columns given by
//
//   U = ( Lambda_{p_1} e_{p_1}, ..., Lambda_{p_M} e_{p_M} )
//
// where e_p is a unit vector in the p-th direction and the
// Lambda_{p_m} are constants.  If p_m appears above, we say the
// direction p_m has been `stabilized.'
//
// We choose which direction to stabilize by comparing diagonal
// entries A_{ii} encountered during the Cholesky decomposition to
// stabilizeThreshold*Lambda_GM, where Lambda_GM is the geometric mean
// of the diagonal entries L computed so far.  Smaller
// stabilizeThreshold means fewer directions will be stabilized.
// Larger stabilizeThreshold means more directions will be stabilized.
//
// Input:
// - A
// - stabilizeThreshold: parameter for deciding which directions to
//   stabilize
// Output:
// - L (overwritten)
// - stabilizeIndices: a list of directions which have been
//   stabilized (overwritten)
// - stabilizeLambdas: a list of corresponding Lambdas (overwritten)
//
void choleskyDecompositionStabilized(Matrix &A, Matrix &L,
                                     vector<Integer> &stabilizeIndices,
                                     vector<Real> &stabilizeLambdas,
                                     const double stabilizeThreshold);

// B := L^{-1} B, where L is lower-triangular and B is a matrix
// pointed to by b
//
// Input:
// - L
// - b, a pointer to the top-left element of B
// - bcols, the number of columns of B
// - ldb, the distance in pointer-space between the first elements of
//   each row of B
// Output:
// - elements of B, which are values pointed to by b are overwritten
//
// (The use of pointers and ldb is to allow using this function for
// submatrices of a larger matrix.)
void lowerTriangularSolve(Matrix &L, Real *b, int bcols, int ldb);

// b := L^{-1} b, where L is lower-triangular
void lowerTriangularSolve(Matrix &L, Vector &b);

// B := L^{-T} B, where L is lower-triangular and B is a matrix
// pointed to by b
//
// Input:
// - L
// - b, a pointer to the top-left element of B
// - bcols, the number of columns of B
// - ldb, the distance in pointer-space between the first elements of
//   each row of B
// Output:
// - elements of B, which are values pointed to by b are overwritten
//
void lowerTriangularTransposeSolve(Matrix &L, Real *b, int bcols, int ldb);

// b := L^{-T} b, where L is lower-triangular
void lowerTriangularTransposeSolve(Matrix &L, Vector &b);

// X := ACholesky^{-1 T} ACholesky^{-1} X = A^{-1} X
// Inputs:
// - ACholesky : dim x dim lower triangular matrix (constant)
// Output:
// - X         : dim x dim matrix (overwritten)
//
void matrixSolveWithCholesky(Matrix &ACholesky, Matrix &X);

#endif  // SDPB_MATRIX_H_
