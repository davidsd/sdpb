#ifndef SDP_BOOTSTRAP_MATRIX_H_
#define SDP_BOOTSTRAP_MATRIX_H_

#include <assert.h>
#include <iostream>
#include <ostream>
#include "omp.h"
#include "types.h"
#include "Vector.h"

using std::ostream;
using std::vector;

const Real CHOLESKY_STABILIZE_THRESHOLD = 1e-8;
const Real BASIC_ROW_THRESHOLD = 1e-4;

class Matrix {
 public:
  int rows;
  int cols;
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

  void setZero() {
    fillVector(elements, 0);
  }

  void addDiagonal(const Real &c) {
    assert(rows == cols);

    for (int i = 0; i < rows; i++)
      elt(i,i) += c;
  }

  void addColumn() {
    elements.resize(elements.size() + rows);
    cols++;
  }

  void setCols(int c) {
    elements.resize(rows*c);
    cols = c;
  }

  void setRowsCols(int r, int c) {
    elements.resize(r*c);
    rows = r;
    cols = c;
  }

  void symmetrize() {
    assert(rows == cols);

    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < r; c++) {
        Real tmp = (elt(r,c) + elt(c,r))/2;
        elt(r,c) = tmp;
        elt(c,r) = tmp;
      }
    }
  }

  void copyFrom(const Matrix &A) {
    assert(rows == A.rows);
    assert(cols == A.cols);

    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] = A.elements[i];
  }

  void operator+=(const Matrix &A) {
    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] += A.elements[i];
  }    

  void operator-=(const Matrix &A) {
    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] -= A.elements[i];
  }

  void operator*=(const Real &c) {
    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] *= c;
  }

  Real maxAbs() const {
    return maxAbsVector(elements);
  }

  void swapCols(int c1, int c2) {
    for (int r = 0; r < rows; r++) {
      Real tmp = elt(r, c1);
      elt(r,c1) = elt(r,c2);
      elt(r,c1) = tmp;
    }
  }

  friend ostream& operator<<(ostream& os, const Matrix& a);
};

ostream& operator<<(ostream& os, const Matrix& a);

// B := A^T
void transpose(const Matrix &A, Matrix &B);

// C := alpha*A*B + beta*C
void matrixScaleMultiplyAdd(Real alpha, Matrix &A, Matrix &B, Real beta, Matrix &C);

// C := A*B
void matrixMultiply(Matrix &A, Matrix &B, Matrix &C);

// B = A^T A
void matrixSquare(Matrix &A, Matrix &B);

// A := L A L^T
void lowerTriangularCongruence(Matrix &A, Matrix &L);

// A := L^{-1} A L^{-T}
void lowerTriangularInverseCongruence(Matrix &A, Matrix &L);

// y := alpha A x + beta y
void vectorScaleMatrixMultiplyAdd(Real alpha, Matrix &A, Vector &x, Real beta, Vector &y);

// y := alpha A^T x + beta y
void vectorScaleMatrixMultiplyTransposeAdd(Real alpha, Matrix &A, Vector &x, Real beta, Vector &y);

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
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
Real minEigenvalue(Matrix &A, Vector &workspace, Vector &eigenvalues);

// Compute an in-place LU decomposition of A, with pivots, suitable
// for use with 'solveWithLUDecomposition'
void LUDecomposition(Matrix &A, vector<Integer> &pivots);

void solveWithLUDecomposition(Matrix &LU, vector<Integer> &pivots, Real *B, int bcols, int ldb);

void solveWithLUDecompositionTranspose(Matrix &LU, vector<Integer> &pivots, Real *B, int bcols, int ldb);

void solveWithLUDecomposition(Matrix &LU, vector<Integer> &pivots, Vector &b);

// L (lower triangular) such that A = L L^T
// Inputs:
// - A : dim x dim symmetric matrix
// - L : dim x dim lower-triangular matrix
void choleskyDecomposition(Matrix &A, Matrix &L);

// L' (lower triangular) such that L' L'^T = L L^T + v v^T. i.e., if L
// is a cholesky decomposition of A, then L' is a cholesky
// decomposition of A + v v^T.  This function dominates the running
// time of this program.
// Inputs: 
// - L : dim x dim lower-triangular matrix 
// - v : pointer to the head of a length-dim vector
// both L and v are modified in place
void choleskyUpdate(Matrix &L, Real *v, int firstNonzeroIndex);

// Let lambdaGeometricMean be the geometric mean of the diagonal
// elements of L.  Whenever a diagonal element is less than
// CHOLESKY_STABILIZE_THRESHOLD * lambdaGeometricMean, update the
// cholesky decomposition L -> L' so that
//
// L' L'^T = L L^T + lambdaGeometricMean^2 u_i u_i^T
//
// where u is a unit vector in the i-th direction.  The index i is
// stored in the vector updateIndices.
void stabilizeCholesky(Matrix &L, Vector &updateVector, vector<int> &updateIndices, Real &lambdaGeometricMean);

void lowerTriangularSolve(Matrix &L, Real *b, int bcols, int ldb);

void lowerTriangularSolve(Matrix &L, Vector &b);

void lowerTriangularTransposeSolve(Matrix &L, Real *b, int bcols, int ldb);

void lowerTriangularTransposeSolve(Matrix &L, Vector &b);

// X := ACholesky^{-1 T} ACholesky^{-1} X = A^{-1} X
// Inputs:
// - ACholesky : dim x dim lower triangular matrix
// - X         : dim x dim matrix
void matrixSolveWithCholesky(Matrix &ACholesky, Matrix &X);

// Return indices of a set of linearly independent rows of M, where
// the basis matrix has eigenvalues whose absolute values exceed
// thresh.
vector<int> linearlyIndependentRowIndices(const Matrix &M);

#endif  // SDP_BOOTSTRAP_MATRIX_H_
