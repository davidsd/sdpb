#include "Matrix.h"

ostream& operator<<(ostream& os, const Matrix& a) {
  os << "{";
  for (int r = 0; r < a.rows; r++) {
    os << "{";
    for (int c = 0; c < a.cols; c++) {
      os << a.elt(r,c);
      if (c < a.cols-1)
        os << ", ";
    }
    os << "}";
    if (r < a.rows-1)
      os << ", ";
  }
  os << "}";
  return os;
}

// B := A^T
void transpose(const Matrix &A, Matrix &B) {
  assert(A.cols == B.rows);
  assert(A.rows == B.cols);

  for (int n = 0; n < A.cols; n++)
    for (int m = 0; m < A.rows; m++)
      B.elt(n,m) = A.elt(m,n);
}

// C := alpha*A*B + beta*C
//
void matrixScaleMultiplyAdd(Real alpha, Matrix &A, Matrix &B, Real beta, Matrix &C) {
  assert(A.cols == B.rows);
  assert(A.rows == C.rows);
  assert(B.cols == C.cols);

  Rgemm("N", "N", A.rows, B.cols, A.cols, alpha,
        &A.elements[0], A.rows,
        &B.elements[0], B.rows,
        beta,
        &C.elements[0], C.rows);
}

// C := A*B
//
void matrixMultiply(Matrix &A, Matrix &B, Matrix &C) {
  matrixScaleMultiplyAdd(1, A, B, 0, C);
}

// B = A^T A
void matrixSquare(Matrix &A, Matrix &B) {
  assert(A.cols == B.rows);
  assert(B.cols == B.rows);
  
  #pragma omp parallel for schedule(dynamic)
  for (int c = 0; c < B.cols; c++) {
    for (int r = 0; r <= c; r++) {
      Real tmp = 0;
      for (int p = 0; p < A.rows; p++)
        tmp += A.elt(p,r) * A.elt(p,c);
      B.elt(r,c) = tmp;
      if (r != c)
        B.elt(c,r) = tmp;
    }
  }
}

// A := L A L^T
void lowerTriangularCongruence(Matrix &A, Matrix &L) {
  int dim = A.rows;
  assert(A.cols == dim);
  assert(L.rows == dim);
  assert(L.cols == dim);

  Rtrmm("Right", "Lower", "Transpose", "NonUnitDiagonal", dim, dim, 1,
        &L.elements[0], dim,
        &A.elements[0], dim);

  Rtrmm("Left", "Lower", "NoTranspose", "NonUnitDiagonal", dim, dim, 1,
        &L.elements[0], dim,
        &A.elements[0], dim);
}

// A := L^{-1} A L^{-T}
void lowerTriangularInverseCongruence(Matrix &A, Matrix &L) {
  int dim = A.rows;
  assert(A.cols == dim);
  assert(L.rows == dim);
  assert(L.cols == dim);

  Rtrsm("Right", "Lower", "Transpose", "NonUnitDiagonal", dim, dim, 1,
        &L.elements[0], dim,
        &A.elements[0], dim);

  Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal", dim, dim, 1,
        &L.elements[0], dim,
        &A.elements[0], dim);
}

// y := alpha A x + beta y
//
void vectorScaleMatrixMultiplyAdd(Real alpha, Matrix &A, Vector &x, Real beta, Vector &y) {
  assert(A.cols == (int)x.size());
  assert(A.rows == (int)y.size());

  #pragma omp parallel for schedule(static)
  for (int p = 0; p < A.rows; p++) {
    Real tmp = 0;
    for (int n = 0; n < A.cols; n++)
      tmp += A.elt(p,n)*x[n];
    y[p] = alpha*tmp + beta*y[p];
  }
}

// y := alpha A^T x + beta y
//
void vectorScaleMatrixMultiplyTransposeAdd(Real alpha, Matrix &A, Vector &x, Real beta, Vector &y) {
  assert(A.cols == (int)y.size());
  assert(A.rows == (int)x.size());

  #pragma omp parallel for schedule(static)
  for (int n = 0; n < A.cols; n++) {
    Real tmp = 0;
    for (int p = 0; p < A.rows; p++)
      tmp += A.elt(p,n)*x[p];
    y[n] = alpha*tmp + beta*y[n];
  }
}

Real frobeniusProductSymmetric(const Matrix &A, const Matrix &B) {
  assert(A.rows == B.rows);
  assert(A.cols == B.cols);
  assert(A.rows == A.cols);

  Real result = 0;
  for (int c = 0; c < A.cols; c++)
    for (int r = 0; r < c ; r++)
      result += A.elt(r,c)*B.elt(r,c);
  result *= 2;

  for (int r = 0; r < A.rows; r++)
    result += A.elt(r,r)*B.elt(r,r);
  
  return result;
}

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric Matrices and
// '.' is the Frobenius product.
//
Real frobeniusProductOfSums(const Matrix &X, const Matrix &dX,
                            const Matrix &Y, const Matrix &dY) {
  Real result = 0;

  for (int c = 0; c < X.cols; c++)
    for (int r = 0; r < c; r++)
      result += (X.elt(r,c) + dX.elt(r,c)) * (Y.elt(r,c) + dY.elt(r,c));
  result *= 2;

  for (int r = 0; r < X.rows; r++)
    result += (X.elt(r,r) + dX.elt(r,r)) * (Y.elt(r,r) + dY.elt(r,r));

  return result;
}

// Eigenvalues of A, via the QR method
// Inputs:
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
//
void matrixEigenvalues(Matrix &A, Vector &workspace, Vector &eigenvalues) {
  assert(A.rows == A.cols);
  assert((int)eigenvalues.size() == A.rows);
  assert((int)workspace.size() == 3*A.rows - 1);

  Integer info;
  Integer workSize = workspace.size();
  Rsyev("NoEigenvectors", "LowerTriangular", A.rows, &A.elements[0], A.rows, &eigenvalues[0], &workspace[0], workSize, &info);
  assert(info == 0);
}

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
//
Real minEigenvalue(Matrix &A, Vector &workspace, Vector &eigenvalues) {
  matrixEigenvalues(A, workspace, eigenvalues);
  return eigenvalues[0];
}

void LUDecomposition(Matrix &A, vector<Integer> &pivots) {
  int dim = A.rows;
  assert(A.cols == dim);

  Integer info;
  Rgetrf(dim, dim, &A.elements[0], dim, &pivots[0], &info);
  assert(info == 0);
}

void solveWithLUDecomposition(Matrix &LU, vector<Integer> &pivots, Real *B, int bcols, int ldb) {
  Integer info;
  Rgetrs("NoTranspose", LU.rows, bcols, &LU.elements[0], LU.rows, &pivots[0], B, ldb, &info);
  assert(info == 0);
}

void solveWithLUDecompositionTranspose(Matrix &LU, vector<Integer> &pivots, Real *B, int bcols, int ldb) {
  Integer info;
  Rgetrs("Transpose", LU.rows, bcols, &LU.elements[0], LU.rows, &pivots[0], B, ldb, &info);
  assert(info == 0);
}

void solveWithLUDecomposition(Matrix &LU, vector<Integer> &pivots, Vector &b) {
  solveWithLUDecomposition(LU, pivots, &b[0], 1, b.size());
}

// L (lower triangular) such that A = L L^T
// Inputs:
// - A : dim x dim symmetric matrix
// - L : dim x dim lower-triangular matrix
//
void choleskyDecomposition(Matrix &A, Matrix &L) {
  int dim = A.rows;
  assert(A.cols == dim);
  assert(L.rows == dim);
  assert(L.cols == dim);

  // Set lower-triangular part of L to cholesky decomposition
  L.copyFrom(A);
  Integer info;
  Rpotrf("Lower", dim, &L.elements[0], dim, &info);
  assert(info == 0);

  // Set the upper-triangular part of the L to zero
  for (int j = 0; j < dim; j++)
    for (int i = 0; i < j; i++)
      L.elements[i + j*dim] = 0;
}

// L' (lower triangular) such that L' L'^T = L L^T + v v^T. i.e., if L
// is a cholesky decomposition of A, then L' is a cholesky
// decomposition of A + v v^T.  This function dominates the running
// time of this program.
// Inputs: 
// - L : dim x dim lower-triangular matrix 
// - v : pointer to the head of a length-dim vector
// both L and v are modified in place
//
void choleskyUpdate(Matrix &L, Real *v, int firstNonzeroIndex) {
  int dim = L.rows;
  Real c, s, x, y;
  for (int r = firstNonzeroIndex; r < dim; r++) {
    x = L.elt(r,r);
    y = *(v+r);
    Rrotg(&x, &y, &c, &s);
    Rrot(dim - r, &L.elements[r*(dim+1)], 1, v+r, 1, c, s);
  }
}

// Let lambdaGeometricMean be the geometric mean of the diagonal
// elements of L.  Whenever a diagonal element is less than
// CHOLESKY_STABILIZE_THRESHOLD * lambdaGeometricMean, update the
// cholesky decomposition L -> L' so that
//
// L' L'^T = L L^T + lambdaGeometricMean^2 u_i u_i^T
//
// where u is a unit vector in the i-th direction.  The index i is
// stored in the vector updateIndices.
//
void stabilizeCholesky(Matrix &L,
                       Vector &updateVector,
                       vector<int> &updateIndices,
                       Real &lambdaGeometricMean) {
  int dim = L.rows;
  assert(L.cols == dim);

  double averageLogLambda = 0;
  for (int i = 0; i < dim; i++) {
    averageLogLambda += log(realToDouble(L.elt(i,i)));
  }
  lambdaGeometricMean = Real(exp(averageLogLambda/dim));
  
  Real lambdaThreshold = CHOLESKY_STABILIZE_THRESHOLD * lambdaGeometricMean;
  for (int i = 0; i < dim; i++) {
    if (L.elt(i,i) < lambdaThreshold) {
      fillVector(updateVector, 0);
      updateVector[i] = lambdaGeometricMean;
      choleskyUpdate(L, &updateVector[0], i);
      updateIndices.push_back(i);
    }
  }
}

void lowerTriangularSolve(Matrix &L, Real *b, int bcols, int ldb) {
  int dim = L.rows;
  assert(L.cols == dim);

  Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal",
        dim, bcols, 1, &L.elements[0], dim, b, ldb);
}

void lowerTriangularSolve(Matrix &L, Vector &b) {
  assert((int) b.size() == L.rows);
  lowerTriangularSolve(L, &b[0], 1, b.size());
}

void lowerTriangularTransposeSolve(Matrix &L, Real *b, int bcols, int ldb) {
  int dim = L.rows;
  assert(L.cols == dim);

  Rtrsm("Left", "Lower", "Transpose", "NonUnitDiagonal",
        dim, bcols, 1, &L.elements[0], dim, b, ldb);
}

void lowerTriangularTransposeSolve(Matrix &L, Vector &b) {
  assert((int) b.size() == L.rows);
  lowerTriangularTransposeSolve(L, &b[0], 1, b.size());
}

// X := ACholesky^{-1 T} ACholesky^{-1} X = A^{-1} X
// Inputs:
// - ACholesky : dim x dim lower triangular matrix
// - X         : dim x dim matrix
//
void matrixSolveWithCholesky(Matrix &ACholesky, Matrix &X) {
  int dim = X.rows;
  assert(X.cols == dim);
  assert(ACholesky.rows == dim);
  assert(ACholesky.cols == dim);

  lowerTriangularSolve(ACholesky, &X.elements[0], X.cols, X.rows);
  lowerTriangularTransposeSolve(ACholesky, &X.elements[0], X.cols, X.rows);
}

// Find the column with the maximum absolute value entry to the right
// of (r,c) (inclusive)
int maxColRightOf(Matrix &A, int r, int c) {
  int  cMax = 0;
  Real aMax = 0;
  for (int k = c; k < A.cols; k++) {
    Real a = abs(A.elt(r,k));
    if (a > aMax) {
      cMax = k;
      aMax = a;
    }
  }
  return cMax;
}

// Return indices of a set of linearly independent rows of M, where
// the basis matrix has eigenvalues whose absolute values exceed
// thresh.
vector<int> linearlyIndependentRowIndices(const Matrix &M) {
  vector<int> rows;
  Matrix A(M);

  int c = 0;
  for (int r = 0; r < A.rows && c < A.cols; r++) {
    A.swapCols(maxColRightOf(A, r, c), c);
    if (abs(A.elt(r, c)) > BASIC_ROW_THRESHOLD) {
      rows.push_back(r);
      // Add a multiple of column c to each column to the right so
      // that the entries to the right of (r,c) are zero.
      for (int k = c+1; k < A.cols; k++) {
        Real x = A.elt(r,k)/A.elt(r,c);
        for (int i = r; i < A.rows; i++)
          A.elt(i,k) -= x*A.elt(i,c);
      }
      c++;
    }
  }
  assert((int)rows.size() == A.cols);

  return rows;
}
