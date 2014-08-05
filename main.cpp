#include <iterator>
#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <assert.h>
#include "omp.h"
#include "types.h"
#include "tinyxml2.h"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/optional.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/timer/timer.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::ostream;
using std::istream;
using std::max;
using std::min;
using std::pair;
using std::ofstream;

using tinyxml2::XMLDocument;
using tinyxml2::XMLElement;

using boost::filesystem::path;
using boost::filesystem::ifstream;
using boost::optional;
using boost::posix_time::second_clock;

using boost::timer::cpu_timer;

class Timers {
public:
  cpu_timer schurComplementCholesky;
  cpu_timer predictorSolution;
  cpu_timer correctorSolution;
  cpu_timer bilinearPairings;
  cpu_timer computeSchurBlocks;
  cpu_timer schurBlocksCholesky;
  cpu_timer schurCholeskyUpdate;
  cpu_timer runSolver;

  Timers() {
    schurComplementCholesky.stop();
    predictorSolution.stop();
    correctorSolution.stop();
    bilinearPairings.stop();
    computeSchurBlocks.stop();
    schurBlocksCholesky.stop();
    schurCholeskyUpdate.stop();
    runSolver.stop();
  }    

  friend ostream& operator<<(ostream& os, const Timers& t);
};

ostream& operator<<(ostream& os, const Timers& t) {
  cout << "Time elapsed:" << endl;
  cout << "  schurComplementCholesky: " << t.schurComplementCholesky.format();
  cout << "  predictorSolution      : " << t.predictorSolution.format();
  cout << "  correctorSolution      : " << t.correctorSolution.format();
  cout << "  bilinearPairings       : " << t.bilinearPairings.format();
  cout << "  computeSchurBlocks     : " << t.computeSchurBlocks.format();
  cout << "  schurBlocksCholesky    : " << t.schurBlocksCholesky.format();
  cout << "  schurCholeskyUpdate    : " << t.schurCholeskyUpdate.format();
  cout << "  runSolver              : " << t.runSolver.format();
  return os;
}

Timers timers;

template <class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  os << "{";
  int last = v.size() - 1;
  for (int i = 0; i < last; i++)
    os << v[i] << ", ";
  if (last >= 0)
    os << v[last];
  os << "}";
  return os;
}

typedef vector<Real> Vector;

Real maxAbsVectorElement(const Vector &v) {
  Real max = 0;
  for (Vector::const_iterator x = v.begin(); x != v.end(); x++)
    if (abs(*x) > max)
      max = abs(*x);
  return max;
}  

void fillVector(Vector &v, const Real &a) {
  std::fill(v.begin(), v.end(), a);
}

// y := alpha*x + beta*y
//
void vectorScaleMultiplyAdd(const Real alpha, const Vector &x, const Real beta, Vector &y) {
  assert(x.size() == y.size());
  
  for (unsigned int i = 0; i < x.size(); i++)
    y[i] = alpha*x[i] + beta*y[i];
}

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

  void setIdentity() {
    assert(rows == cols);

    setZero();
    addDiagonal(1);
  }

  void symmetrize() {
    assert(rows == cols);

    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < r; c++) {
        elt(r,c) /= 2;
        elt(r,c) += elt(c,r)/2;
        elt(c,r) = elt(r,c);
        //set(r, c, tmp);
        //set(c, r, tmp);
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

  Real maxAbsElement() const {
    return maxAbsVectorElement(elements);
  }

  friend ostream& operator<<(ostream& os, const Matrix& a);
};

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

Real dotProduct(const Vector &u, const Vector v) {
  Real result = 0;
  for (unsigned int i = 0; i < u.size(); i++)
    result += u[i]*v[i];
  return result;
}

// y := alpha A x + beta y
//
void vectorScaleMatrixMultiplyAdd(Real alpha, Matrix &A, Vector &x, Real beta, Vector &y) {
  assert(A.cols == (int)x.size());
  assert(A.rows == (int)y.size());

  Rgemv("NoTranspose",
        A.rows, A.cols, alpha,
        &A.elements[0], A.rows,
        &x[0], 1,
        beta,
        &y[0], 1);
}

// y := alpha A^T x + beta y
//
void vectorScaleMatrixMultiplyTransposeAdd(Real alpha, Matrix &A, Vector &x, Real beta, Vector &y) {
  assert(A.cols == (int)y.size());
  assert(A.rows == (int)x.size());

  Rgemv("Transpose",
        A.rows, A.cols, alpha,
        &A.elements[0], A.rows,
        &x[0], 1,
        beta,
        &y[0], 1);
}

void lowerTriangularMatrixTimesVector(Matrix &A, Vector &v) {
  int dim = A.rows;
  assert(A.cols == dim);
  assert((int)v.size() == dim);
  Rtrmv("Lower", "NoTranspose", "NotUnitDiagonal", dim, &A.elements[0], dim, &v[0], 1);
}

void lowerTriangularMatrixTransposeTimesVector(Matrix &A, Vector &v) {
  int dim = A.rows;
  assert(A.cols == dim);
  assert((int)v.size() == dim);
  Rtrmv("Lower", "Transpose", "NotUnitDiagonal", dim, &A.elements[0], dim, &v[0], 1);
}

Real frobeniusProduct(const Matrix &A, const Matrix &B) {
  assert(A.rows == B.rows);
  assert(A.cols == B.cols);
  return dotProduct(A.elements, B.elements);
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

void LUDecomposition(Matrix &A, vector<mpackint> &ipiv) {
  int dim = A.rows;
  assert(A.cols == dim);

  mpackint info;
  Rgetrf(dim, dim, &A.elements[0], dim, &ipiv[0], &info);
  assert(info == 0);
}

void solveWithLUDecomposition(Matrix &LU, vector<mpackint> &ipiv, Real *B, int bcols, int ldb) {
  mpackint info;
  Rgetrs("NoTranspose", LU.rows, bcols, &LU.elements[0], LU.rows, &ipiv[0], B, ldb, &info);
  assert(info == 0);
}

void solveWithLUDecompositionTranspose(Matrix &LU, vector<mpackint> &ipiv, Real *B, int bcols, int ldb) {
  mpackint info;
  Rgetrs("Transpose", LU.rows, bcols, &LU.elements[0], LU.rows, &ipiv[0], B, ldb, &info);
  assert(info == 0);
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

  if (dim == 1) {
    L.elt(0, 0) = sqrt(A.elt(0, 0));
    return;
  }

  // Set lower-triangular part of L to cholesky decomposition
  L.copyFrom(A);
  mpackint info;
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
void choleskyUpdate(Matrix &L, Real *v) {
  int dim = L.rows;
  Real c, s, x, y;
  for (int r = 0; r < dim; r++) {
    x = L.elt(r,r);
    y = *(v+r);
    Rrotg(&x, &y, &c, &s);

    Real *dx = &L.elements[r*(dim+1)];
    Real *dy = v+r;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < dim - r; i++) {
      const Real tmp = c*dx[i] + s*dy[i];
      dy[i] = c*dy[i] - s*dx[i];
      dx[i] = tmp;
    }
  }
}

// L' (lower triangular) such that L' L'^T = L L^T + V V^T. i.e., if L
// is a cholesky decomposition of A, then L' is a cholesky
// decomposition of A + V V^T.  This is more efficient than directly
// computing the cholesky decomposition of A + V V^T if V has a small
// number of columns.
// Inputs: 
// - L : dim x dim lower-triangular matrix 
// - V : dim x n matrix
// both L and V are modified in place
//
void choleskyUpdate(Matrix &L, Matrix &V) {
  assert(L.rows == V.rows);
  for (int c = 0; c < V.cols; c++)
    choleskyUpdate(L, &V.elements[c*V.rows]);
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

// b := ACholesky^{-1 T} ACholesky^{-1} b = A^{-1} b
//
// Inputs:
// - ACholesky : dim x dim lower triangular matrix, the Cholesky decomposition of a matrix A
// - b         : vector of length dim (output)
//
void vectorSolveWithCholesky(Matrix &ACholesky, Vector &b) {
  assert((int) b.size() == ACholesky.rows);

  lowerTriangularSolve(ACholesky, &b[0], 1, b.size());
  lowerTriangularTransposeSolve(ACholesky, &b[0], 1, b.size());
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

// result = b'^T a b', where b' = b \otimes 1
// Inputs:
// - a      : l*m x l*m symmetric matrix
// - b      : l   x n   matrix
// - work   : l*m x n*m matrix
// - result : n*m x n*m symmetric matrix
//
void tensorMatrixCongruenceTranspose(const Matrix &a,
                                     const Matrix &b,
                                     Matrix &work,
                                     Matrix &result) {
  int m = a.rows / b.rows;

  assert(result.rows == b.cols * m);
  assert(result.cols == b.cols * m);

  assert(work.rows == a.rows);
  assert(work.cols == result.cols);

  // work = a b'
  for (int c = 0; c < work.cols; c++) {
    int bCol       = c % b.cols;
    int aColOffset = (c / b.cols) * b.rows;

    for (int r = 0; r < work.rows; r++) {

      Real tmp = 0;
      for (int k = 0; k < b.rows; k++) {
        tmp += a.elt(r, aColOffset + k) * b.elt(k, bCol);
      }

      work.elt(r, c) = tmp;
    }
  }

  // result = b'^T work
  for (int c = 0; c < result.cols; c++) {

    // since result is symmetric, only compute its upper triangle
    for (int r = 0; r <= c; r++) {
      int bCol          = r % b.cols;
      int workRowOffset = (r / b.cols) * b.rows;

      Real tmp = 0;
      for (int k = 0; k < b.rows; k++) {
        tmp += b.elt(k, bCol) * work.elt(workRowOffset + k, c);
      }

      result.elt(r, c) = tmp;

      // lower triangle is the same as upper triangle
      if (c != r) {
        result.elt(c, r) = tmp;
      }
    }
  }
}

// result = B'^T A^{-1} B', where B' = B \otimes 1
// Inputs:
// - L      : l*m x l*m cholesky decomposition of A
// - B      : l   x n   matrix
// - Work   : l*m x n*m matrix
// - Result : n*m x n*m symmetric matrix
//
void tensorMatrixInvCongruenceTransposeWithCholesky(const Matrix &L,
                                                    const Matrix &B,
                                                    Matrix &Work,
                                                    Matrix &Result) {
  // X = L^{-1} (B \otimes 1);
  for (int cw = 0; cw < Work.cols; cw++) {
    int mc  = cw / B.cols;

    for (int rw = mc*B.rows; rw < Work.rows; rw++) {
      int mr = rw / B.cols;

      Real tmp = (mr != mc) ? 0 : B.elt(rw % B.rows, cw % B.cols);
      for (int cl = 0; cl < rw; cl++)
        tmp -= L.elt(rw, cl)*Work.elt(cl, cw);

      Work.elt(rw, cw) = tmp/L.elt(rw, rw);
    }
  }

  // // Work = LInv (B \otimes 1)
  // for (int cw = 0; cw < Work.cols; cw++) {
  //   int m  = cw / B.cols;
  //   int cb = cw % B.cols;

  //   for (int rw = m*B.rows; rw < Work.rows; rw++) {
  //     Real tmp = 0;
  //     for (int cl = m*B.rows; cl < min(rw+1, (m+1)*B.rows); cl++)
  //       tmp += LInv.elt(rw, cl)*B.elt(cl % B.rows, cb);
  //     Work.elt(rw, cw) = tmp;
  //   }
  // }

  // Result = Work^T Work
  for (int cr = 0; cr < Result.cols; cr++) {
    int mc = cr / B.cols;

    for (int rr = 0; rr <= cr; rr++) {
      int mr = rr / B.cols;

      Real tmp = 0;
      for (int rw = max(mr, mc)*B.rows; rw < Work.rows; rw++)
        tmp += Work.elt(rw, cr)*Work.elt(rw, rr);

      Result.elt(rr, cr) = tmp;
      if (rr != cr)
        Result.elt(cr, rr) = tmp;
    }
  }
}

void testTensorCongruence() {
  Matrix a(4,4);
  Matrix L(a);
  Matrix b(2,3);
  Matrix result(6,6);
  Matrix work(4,6);
  a.setIdentity();
  choleskyDecomposition(a,L);
  b.elt(0,0) =2;
  b.elt(1,0) =3;
  b.elt(0,1) =4;
  b.elt(1,1) =5;
  b.elt(0,2) =6;
  b.elt(1,2) =7;

  tensorMatrixCongruenceTranspose(a, b, work, result);
  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  cout << "work = " << work << endl;
  cout << "result = " << result << endl;

  tensorMatrixInvCongruenceTransposeWithCholesky(L, b, work, result);
  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  cout << "work = " << work << endl;
  cout << "result = " << result << endl;

}

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

  void setIdentity() {
    setZero();
    addDiagonal(1);
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

  Real maxAbsElement() const {
    Real max = 0;
    for (vector<Matrix>::const_iterator b = blocks.begin(); b != blocks.end(); b++) {
      const Real tmp = b->maxAbsElement();
      if (tmp > max)
        max = tmp;
    }
    return max;
  }

  friend ostream& operator<<(ostream& os, const BlockDiagonalMatrix& A);

};

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
  // TODO: Parallelize this loop
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    result += frobeniusProductSymmetric(A.blocks[b], B.blocks[b]);
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
  // TODO: Parallelize this loop
  for (unsigned int b = 0; b < X.blocks.size(); b++)
    result += frobeniusProductOfSums(X.blocks[b], dX.blocks[b], Y.blocks[b], dY.blocks[b]);
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

void choleskyDecomposition(BlockDiagonalMatrix &A,
                           BlockDiagonalMatrix &L) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    choleskyDecomposition(A.blocks[b], L.blocks[b]);
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

class IndexTuple {
public:
  int p;
  int r;
  int s;
  int k;
  IndexTuple(int p, int r, int s, int k): p(p), r(r), s(s), k(k) {}
  IndexTuple() {}
};

class SDP {
public:
  vector<Matrix> bilinearBases;
  Matrix polMatrixValues;
  Vector affineConstants;
  Vector objective;
  vector<int> dimensions;
  vector<int> degrees;
  vector<vector<int> > blocks;
  vector<vector<IndexTuple> > constraintIndices;

  int numConstraints() const {
    return polMatrixValues.rows;
  }

  vector<int> psdMatrixBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      for (vector<int>::const_iterator b = blocks[j].begin(); b != blocks[j].end(); b++)
        dims.push_back(bilinearBases[*b].rows * dimensions[j]);
    return dims;
  }

  vector<int> bilinearPairingBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      for (vector<int>::const_iterator b = blocks[j].begin(); b != blocks[j].end(); b++)
        dims.push_back(bilinearBases[*b].cols * dimensions[j]);
    return dims;
  }

  vector<int> schurBlockDims() const {
    vector<int> dims;
    for (unsigned int j = 0; j < dimensions.size(); j++)
      dims.push_back(constraintIndices[j].size());
    return dims;
  }

  void initializeConstraintIndices() {
    int p = 0;
    for (unsigned int j = 0; j < dimensions.size(); j++) {
      constraintIndices.push_back(vector<IndexTuple>(0));

      for (int s = 0; s < dimensions[j]; s++) {
        for (int r = 0; r <= s; r++) {
          for (int k = 0; k <= degrees[j]; k++) {
            constraintIndices[j].push_back(IndexTuple(p, r, s, k));
            p++;
          }
        }
      }
    }
    assert(p == numConstraints());
  }

  friend ostream& operator<<(ostream& os, const SDP& sdp);
};

ostream& operator<<(ostream& os, const SDP& sdp) {
  os << "SDP(bilinearBases = " << sdp.bilinearBases
     << ", polMatrixValues = " << sdp.polMatrixValues
     << ", affineConstants = " << sdp.affineConstants
     << ", objective = " << sdp.objective
     << ", dimensions = " << sdp.dimensions
     << ", degrees = " << sdp.degrees
     << ", blocks = " << sdp.blocks
     << ")";

  return os;
}

template <class T>
vector<T> parseMany(const char *name, T(*parse)(XMLElement *), XMLElement *elt) {
  XMLElement *e;
  vector<T> v;
  for (e = elt->FirstChildElement(name);
       e != NULL;
       e = e->NextSiblingElement(name)) {
    v.push_back(parse(e));
  }
  return v;
}

Real parseReal(XMLElement *rXml) {
  return Real(rXml->GetText());
}

int parseInt(XMLElement *iXml) {
  return atoi(iXml->GetText());
}

Vector parseVector(XMLElement *vecXml) {
  return parseMany("elt", parseReal, vecXml);
}

Matrix parseMatrix(XMLElement *matXml) {
  Matrix m;
  m.rows     = parseInt(matXml->FirstChildElement("rows"));
  m.cols     = parseInt(matXml->FirstChildElement("cols"));
  m.elements = parseVector(matXml->FirstChildElement("elements"));
  return m;
}

class SampledMatrixPolynomial {
public:
  int dim;
  int degree;
  Matrix samplesMatrix;
  vector<Matrix> bilinearBases;
};

SampledMatrixPolynomial parseSampledMatrixPolynomial(XMLElement *smpXml) {
  SampledMatrixPolynomial s;
  s.dim    = parseInt(smpXml->FirstChildElement("dim"));
  s.degree = parseInt(smpXml->FirstChildElement("degree"));
  s.samplesMatrix = parseMatrix(smpXml->FirstChildElement("samplesMatrix"));
  s.bilinearBases = parseMany("bilinearBasisMatrix", parseMatrix, smpXml->FirstChildElement("bilinearBases"));
  return s;
}

SDP bootstrapSDP(const Vector &objective,
                 const Vector &normalization,
                 const vector<SampledMatrixPolynomial> &sampledMatrixPols) {
  SDP sdp;
  sdp.objective = objective;

  int numConstraints = 0;
  for (vector<SampledMatrixPolynomial>::const_iterator s = sampledMatrixPols.begin();
       s != sampledMatrixPols.end();
       s++) {
    int dimension = s->dim;
    int degree    = s->degree;

    sdp.dimensions.push_back(dimension);
    sdp.degrees.push_back(degree);
    numConstraints += (degree+1)*dimension*(dimension+1)/2;
  }

  // For the normalization constraint
  // sdp.dimensions.push_back(1);
  // sdp.degrees.push_back(0);
  // numConstraints += 1;

  sdp.polMatrixValues = Matrix(numConstraints, sdp.objective.size());
  // sdp.affineConstants = Vector(numConstraints, 0);
  // THIS IS A HACK!!!
  sdp.affineConstants = Vector(numConstraints, 1);

  // normalization constraint
  // sdp.affineConstants[numConstraints-1] = 1;

  int p = 0;
  for (vector<SampledMatrixPolynomial>::const_iterator s = sampledMatrixPols.begin();
       s != sampledMatrixPols.end();
       s++) {

    vector<int> blocks;
    for (vector<Matrix>::const_iterator b = s->bilinearBases.begin();
         b != s->bilinearBases.end();
         b++) {
      assert(b->cols == s->degree + 1);
      blocks.push_back(sdp.bilinearBases.size());
      sdp.bilinearBases.push_back(*b);
    }
    sdp.blocks.push_back(blocks);

    for (int k = 0; k < s->samplesMatrix.rows; k++, p++)
      for (int n = 0; n < s->samplesMatrix.cols; n++)
        sdp.polMatrixValues.elt(p, n) = -(s->samplesMatrix.elt(k, n));
  }
  // assert(p == numConstraints - 1);
  assert(p == numConstraints);

  // normalization constraint
  // for (unsigned int n = 0; n < sdp.objective.size(); n++)
  //   sdp.polMatrixValues.elt(p, n) = normalization[n];
  // sdp.blocks.push_back(vector<int>());

  //sdp.rescale();
  // sdp.addSlack();
  sdp.initializeConstraintIndices();
  return sdp;
}


SDP parseBootstrapSDP(XMLElement *sdpXml) {
  return bootstrapSDP(parseVector(sdpXml->FirstChildElement("objective")),
                      parseVector(sdpXml->FirstChildElement("normalization")),
                      parseMany("sampledMatrixPolynomial",
                                parseSampledMatrixPolynomial,
                                sdpXml->FirstChildElement("sampledPositiveMatrices")));
}

SDP readBootstrapSDP(const path sdpFile) {
  XMLDocument doc;
  doc.LoadFile(sdpFile.c_str());
  return parseBootstrapSDP(doc.FirstChildElement("sdp"));
}

class SDPSolverParameters {
public:
  int maxIterations;
  Real dualityGapThreshold;          // = epsilonStar
  Real primalFeasibilityThreshold;   // = epsilonBar
  Real dualFeasibilityThreshold;     // = epsilonBar
  Real initialMatrixScale;           // = lambdaStar
  Real feasibleCenteringParameter;   // = betaStar
  Real infeasibleCenteringParameter; // = betaBar
  Real stepLengthReduction;          // = gammaStar
  int precision;
  int maxThreads;

  SDPSolverParameters():
    maxIterations(130),
    dualityGapThreshold("1e-100"),
    primalFeasibilityThreshold("1e-30"),
    dualFeasibilityThreshold("1e-30"),
    initialMatrixScale("1e20"),
    feasibleCenteringParameter("0.1"),
    infeasibleCenteringParameter("0.3"),
    stepLengthReduction("0.7"),
    precision(500),
    maxThreads(1) {}

  SDPSolverParameters(const path &paramFile) {
    ifstream ifs(paramFile);
    ifs >> maxIterations;                ifs.ignore(256, '\n');
    ifs >> dualityGapThreshold;          ifs.ignore(256, '\n');
    ifs >> primalFeasibilityThreshold;   ifs.ignore(256, '\n');
    ifs >> dualFeasibilityThreshold;     ifs.ignore(256, '\n');
    ifs >> initialMatrixScale;           ifs.ignore(256, '\n');
    ifs >> feasibleCenteringParameter;   ifs.ignore(256, '\n');
    ifs >> infeasibleCenteringParameter; ifs.ignore(256, '\n');
    ifs >> stepLengthReduction;          ifs.ignore(256, '\n');
    ifs >> precision;                    ifs.ignore(256, '\n');
    ifs >> maxThreads;                   ifs.ignore(256, '\n');
  }

  friend ostream& operator<<(ostream& os, const SDPSolverParameters& p);

};

ostream& operator<<(ostream& os, const SDPSolverParameters& p) {
  os << "maxIterations                = " << p.maxIterations                << endl;
  os << "dualityGapThreshold          = " << p.dualityGapThreshold          << endl;
  os << "primalFeasibilityThreshold   = " << p.primalFeasibilityThreshold   << endl;
  os << "dualFeasibilityThreshold     = " << p.dualFeasibilityThreshold     << endl;
  os << "initialMatrixScale           = " << p.initialMatrixScale           << endl;
  os << "feasibleCenteringParameter   = " << p.feasibleCenteringParameter   << endl;
  os << "infeasibleCenteringParameter = " << p.infeasibleCenteringParameter << endl;
  os << "stepLengthReduction          = " << p.stepLengthReduction          << endl;
  os << "precision(actual)            = " << p.precision << "(" << mpf_get_default_prec() << ")" << endl;
  os << "maxThreads                   = " << p.maxThreads                   << endl;
  return os;
}

class SDPSolverStatus {
public:
  Real primalObjective;
  Real dualObjective;
  Real primalError;
  Real dualError;

  Real dualityGap() const {
    return abs(primalObjective - dualObjective) /
      max(Real(abs(primalObjective) + abs(dualObjective)), Real(1));
  }

  bool isPrimalFeasible(const SDPSolverParameters &p) {
    return primalError < p.primalFeasibilityThreshold;
  }

  bool isDualFeasible(const SDPSolverParameters &p) {
    return dualError < p.dualFeasibilityThreshold;
  }

  bool isOptimal(const SDPSolverParameters &p) {
    return dualityGap() < p.dualityGapThreshold;
  }

  friend ostream& operator<<(ostream& os, const SDPSolverStatus& s);
  
};

ostream& operator<<(ostream& os, const SDPSolverStatus& s) {
  os << "primalObjective = " << s.primalObjective << endl;
  os << "dualObjective   = " << s.dualObjective << endl;
  os << "dualityGap      = " << s.dualityGap() << endl;
  os << "primalError     = " << s.primalError << endl;
  os << "dualError       = " << s.dualError << endl;
  return os;
}

class SDPSolver {
public:
  SDP sdp;
  SDPSolverStatus status;

  // current point
  Vector x;
  BlockDiagonalMatrix X;
  BlockDiagonalMatrix Y;

  // search direction
  Vector dx;
  BlockDiagonalMatrix dX;
  BlockDiagonalMatrix dY;

  // discrepancies in dual and primal equality constraints
  Vector dualResidues;
  Vector dualResiduesTilde;
  BlockDiagonalMatrix PrimalResidues;

  // For free variable elimination
  Matrix E;
  Vector g;

  // intermediate computations
  BlockDiagonalMatrix XCholesky;
  BlockDiagonalMatrix YCholesky;
  BlockDiagonalMatrix Z;
  BlockDiagonalMatrix R;
  BlockDiagonalMatrix BilinearPairingsXInv;
  BlockDiagonalMatrix BilinearPairingsY;
  BlockDiagonalMatrix SchurBlocks;
  BlockDiagonalMatrix SchurBlocksCholesky;
  Matrix SchurUpdateLowRank;
  Matrix Q;
  Matrix M;
  Vector basicKernelCoords;
  Matrix BasicKernelSpan;

  // additional workspace variables
  BlockDiagonalMatrix StepMatrixWorkspace;
  vector<Matrix> bilinearPairingsWorkspace;
  vector<Vector> eigenvaluesWorkspace;
  vector<Vector> QRWorkspace;

  SDPSolver(const SDP &sdp):
    sdp(sdp),
    x(sdp.numConstraints(), 0),
    X(sdp.psdMatrixBlockDims()),
    Y(X),
    dx(x),
    dX(X),
    dY(Y),
    dualResidues(x),
    dualResiduesTilde(sdp.numConstraints() - sdp.objective.size()),
    PrimalResidues(X),
    E(sdp.numConstraints() - sdp.objective.size(), sdp.objective.size()),
    g(sdp.objective.size()),
    XCholesky(X),
    YCholesky(X),
    Z(X),
    R(X),
    BilinearPairingsXInv(sdp.bilinearPairingBlockDims()),
    BilinearPairingsY(BilinearPairingsXInv),
    SchurBlocks(sdp.schurBlockDims()),
    SchurBlocksCholesky(SchurBlocks),
    SchurUpdateLowRank(sdp.polMatrixValues),
    Q(sdp.polMatrixValues.cols, sdp.polMatrixValues.cols),
    M(Q),
    basicKernelCoords(Q.rows),
    BasicKernelSpan(sdp.polMatrixValues),
    StepMatrixWorkspace(X)
  {
    // initialize bilinearPairingsWorkspace, eigenvaluesWorkspace, QRWorkspace 
    for (unsigned int b = 0; b < sdp.bilinearBases.size(); b++) {
      bilinearPairingsWorkspace.push_back(Matrix(X.blocks[b].rows, BilinearPairingsXInv.blocks[b].cols));
      eigenvaluesWorkspace.push_back(Vector(X.blocks[b].rows));
      QRWorkspace.push_back(Vector(3*X.blocks[b].rows - 1));
    }

    // Computations needed for free variable elimination
    Matrix DBLU(sdp.objective.size(), sdp.objective.size());
    vector<mpackint> DBLUipiv(sdp.objective.size());

    // LU Decomposition of D_B
    for (int n = 0; n < DBLU.cols; n++)
      for (int m = 0; m < DBLU.rows; m++)
        DBLU.elt(m,n) = sdp.polMatrixValues.elt(m,n);
    LUDecomposition(DBLU, DBLUipiv);

    // Compute E = - D_N D_B^{-1}
    int nonBasicStart = sdp.objective.size();
    // ET = -D_N^T
    Matrix ET(E.cols, E.rows);
    for (int p = 0; p < ET.cols; p++)
      for (int n = 0; n < ET.rows; n++)
        ET.elt(n, p) = -sdp.polMatrixValues.elt(p + nonBasicStart, n);
    // ET = D_B^{-1 T} ET = -D_B^{-1 T} D_N^T
    solveWithLUDecompositionTranspose(DBLU, DBLUipiv, &ET.elements[0], ET.cols, ET.rows);
    // E = ET^T
    transpose(ET, E);

    // g = D_B^{-T} f
    for (unsigned int n = 0; n < g.size(); n++)
      g[n] = sdp.objective[n];
    solveWithLUDecompositionTranspose(DBLU, DBLUipiv, &g[0], 1, g.size());

    // BasicKernelSpan = ( -1 \\ E)
    BasicKernelSpan.setZero();
    for (int c = 0; c < E.cols; c++)
      for (int r = 0; r < E.rows; r++)
        BasicKernelSpan.elt(nonBasicStart + r, c) = E.elt(r, c);
    for (int c = 0; c < E.cols; c++)
      BasicKernelSpan.elt(c, c) = -1;
  }

  void initialize(const SDPSolverParameters &parameters);
  SDPSolverStatus run(const SDPSolverParameters &parameters, const path outFile, const path checkpointFile);
  void computeSearchDirectionWithRMatrix(const BlockDiagonalMatrix &R, const bool isPrimalFeasible);

};

void computeBilinearPairings(const BlockDiagonalMatrix &A,
                             const vector<Matrix> &bilinearBases,
                             vector<Matrix> &workspace,
                             BlockDiagonalMatrix &result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < bilinearBases.size(); b++)
    tensorMatrixCongruenceTranspose(A.blocks[b], bilinearBases[b], workspace[b], result.blocks[b]);
}

void computeInvBilinearPairingsWithCholesky(const BlockDiagonalMatrix &L,
                                            const vector<Matrix> &bilinearBases,
                                            vector<Matrix> &workspace,
                                            BlockDiagonalMatrix &result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < bilinearBases.size(); b++)
    tensorMatrixInvCongruenceTransposeWithCholesky(L.blocks[b], bilinearBases[b], workspace[b], result.blocks[b]);
}

// result = V D V^T, where D=diag(d) is a diagonal matrix
// Inputs:
// - d        : pointer to beginning of a length-V.cols vector
// - V        : V.rows x V.cols Matrix
// - blockRow : integer < k
// - blockCol : integer < k
// - result   : (k*V.rows) x (k*V.rows) square Matrix
//
void diagonalCongruence(Real const *d,
                        const Matrix &V,
                        const int blockRow,
                        const int blockCol,
                        Matrix &result) {

  for (int p = 0; p < V.rows; p++) {
    for (int q = 0; q <= p; q++) {
      Real tmp = 0;

      for (int n = 0; n < V.cols; n++)
        tmp += *(d+n) * V.elt(p, n)*V.elt(q, n);
      
      result.elt(blockRow*V.rows + p, blockCol*V.rows + q) = tmp;
      if (p != q)
        result.elt(blockRow*V.rows + q, blockCol*V.rows + p) = tmp;
    }
  }
}

// v^T A' v, where A' is the (blockRow,blockCol)-th dim x dim block
// inside A.
//
// Input:
// - v        : pointer to the beginning of a vector of length dim
// - dim      : length of the vector v
// - A        : (k*dim) x (k*dim) matrix, where k > blockRow, blockCol
// - blockRow : integer labeling block of A
// - blockCol : integer labeling block of A
//
Real bilinearBlockPairing(const Real *v,
                          const int dim,
                          const Matrix &A,
                          const int blockRow,
                          const int blockCol) {
  Real result = 0;

  for (int r = 0; r < dim; r++) {
    Real tmp = 0;

    for (int c = 0; c < dim; c++)
      tmp += *(v+c) * A.elt(blockRow*dim + r, blockCol*dim + c);
    result += *(v+r) * tmp;
  }
  return result;
}

void computeSchurBlocks(const SDP &sdp,
                        const BlockDiagonalMatrix &BilinearPairingsXInv,
                        const BlockDiagonalMatrix &BilinearPairingsY,
                        BlockDiagonalMatrix &SchurBlocks) {

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] + 1;

    for (unsigned int u1 = 0; u1 < sdp.constraintIndices[j].size(); u1++) {
      const int ej_r1 = sdp.constraintIndices[j][u1].r * ej;
      const int ej_s1 = sdp.constraintIndices[j][u1].s * ej;
      const int k1    = sdp.constraintIndices[j][u1].k;

      for (unsigned int u2 = 0; u2 < sdp.constraintIndices[j].size(); u2++) {
        const int ej_r2 = sdp.constraintIndices[j][u2].r * ej;
        const int ej_s2 = sdp.constraintIndices[j][u2].s * ej;
        const int k2    = sdp.constraintIndices[j][u2].k;

        Real tmp = 0;
        for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
          tmp += (BilinearPairingsXInv.blocks[*b].elt(ej_s1 + k1, ej_r2 + k2) * BilinearPairingsY.blocks[*b].elt(ej_s2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_r1 + k1, ej_r2 + k2) * BilinearPairingsY.blocks[*b].elt(ej_s2 + k2, ej_s1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_s1 + k1, ej_s2 + k2) * BilinearPairingsY.blocks[*b].elt(ej_r2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_r1 + k1, ej_s2 + k2) * BilinearPairingsY.blocks[*b].elt(ej_r2 + k2, ej_s1 + k1))/4;
        }
        SchurBlocks.blocks[j].elt(u1, u2) = tmp;
        if (u2 != u1)
          SchurBlocks.blocks[j].elt(u2, u1) = tmp;
      }
    }
  }
}

// x_B = g + E^T x_N
void basicCompletion(const Vector &g, const Matrix &E, Vector &x) {
  int nonBasicStart = E.cols;
  assert((int)x.size() == E.cols + E.rows);
  
  for (int n = 0; n < nonBasicStart; n++) {
    x[n] = g[n];
    for (int p = 0; p < E.rows; p++)
      x[n] += E.elt(p, n) * x[nonBasicStart + p];
  }
}

// xTilde_N = x_N + E x_B
void nonBasicShift(const Matrix &E, const Vector &x, Vector &xTilde) {
  int nonBasicStart = x.size() - xTilde.size();
  assert(E.rows == (int)xTilde.size());
  assert(E.cols == nonBasicStart);
  
  for (unsigned int p = 0; p < xTilde.size(); p++) {
    xTilde[p] = x[p + nonBasicStart];
    for (int n = 0; n < E.cols; n++)
      xTilde[p] += E.elt(p,n) * x[n];
  }
}

void computeDualResidues(const SDP &sdp,
                         const BlockDiagonalMatrix &Y,
                         const BlockDiagonalMatrix &BilinearPairingsY,
                         Vector &dualResidues) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] +1;

    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t++) {
      const int p    = t->p;
      const int ej_r = t->r * ej;
      const int ej_s = t->s * ej;
      const int k    = t->k;

      dualResidues[p] = 0;
      for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
        dualResidues[p] -= BilinearPairingsY.blocks[*b].elt(ej_r+k, ej_s+k);
        dualResidues[p] -= BilinearPairingsY.blocks[*b].elt(ej_s+k, ej_r+k);
      }
      dualResidues[p] /= 2;
      dualResidues[p] += sdp.affineConstants[p];
    }
  }
}

void constraintMatrixWeightedSum(const SDP &sdp, const Vector x, BlockDiagonalMatrix &result)  {
  // TODO: parallelize this loop
  int p = 0;
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int dj = sdp.degrees[j];

    for (int s = 0; s < sdp.dimensions[j]; s++) {
      for (int r = 0; r <= s; r++) {
        for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++)
          diagonalCongruence(&x[p], sdp.bilinearBases[*b], r, s, result.blocks[*b]);
        p += dj + 1;
      }
    }
  }
  assert(p == (int)x.size());

  result.symmetrize();
}

void computeSchurRHS(const SDP &sdp,
                     Vector &dualResidues,
                     BlockDiagonalMatrix &Z, 
                     Vector &r) {

  for (unsigned int p = 0; p < r.size(); p++)
    r[p] = -dualResidues[p];

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t++) {
      for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {

        const int delta = sdp.bilinearBases[*b].rows;
        // Pointer to the k-th column of sdp.bilinearBases[*b]
        const Real *q = &sdp.bilinearBases[*b].elements[(t->k) * delta];

        r[t->p] -= bilinearBlockPairing(q, delta, Z.blocks[*b], t->r, t->s);
      }      
    }
  }
}

void SDPSolver::initialize(const SDPSolverParameters &parameters) {
  fillVector(x, 0);
  X.setZero();
  X.addDiagonal(parameters.initialMatrixScale);
  Y.setZero();
  Y.addDiagonal(parameters.initialMatrixScale);
}

// PrimalResidues = sum_p F_p x_p - X - F_0
//
void computePrimalResidues(const SDP &sdp,
                           const Vector x,
                           const BlockDiagonalMatrix &X,
                           BlockDiagonalMatrix &PrimalResidues) {
  constraintMatrixWeightedSum(sdp, x, PrimalResidues);
  PrimalResidues -= X;
  // PrimalResidues.addDiagonalPart(sdp.objective, -1);
}

Real primalObjective(const SDP &sdp, const Vector &x) {
  return dotProduct(sdp.affineConstants, x);
}

Real dualObjective(const SDP &sdp, const Vector &g, const Vector &dualResidues) {
  Real tmp = 0;
  for (unsigned int i = 0; i < g.size(); i++)
    tmp += g[i]*dualResidues[i];
  return tmp;
}

// Implements SDPA's DirectionParameter::MehrotraPredictor
Real predictorCenteringParameter(const SDPSolverParameters &parameters, 
                                 const bool reductionSwitch,
                                 const bool isPrimalDualFeasible) {
  if (isPrimalDualFeasible)
    return 0;
  else if (reductionSwitch)
    return parameters.infeasibleCenteringParameter;
  else
    return 2;
}

// Implements SDPA's DirectionParameter::MehrotraCorrector
Real correctorCenteringParameter(const SDPSolverParameters &parameters,
                                 const BlockDiagonalMatrix &X,
                                 const BlockDiagonalMatrix &dX,
                                 const BlockDiagonalMatrix &Y,
                                 const BlockDiagonalMatrix &dY,
                                 const Real &mu,
                                 const bool isPrimalDualFeasible) {

  Real r = frobeniusProductOfSums(X, dX, Y, dY) / (mu * X.dim);
  Real beta = r < 1 ? r*r : r;

  if (isPrimalDualFeasible)
    return min(max(parameters.feasibleCenteringParameter, beta), Real(1));
  else
    return max(parameters.infeasibleCenteringParameter, beta);
}

// R = beta mu I - X Y
//
void computePredictorRMatrix(const Real &beta,
                             const Real &mu,
                             BlockDiagonalMatrix &X,
                             BlockDiagonalMatrix &Y,
                             BlockDiagonalMatrix &R) {
  blockDiagonalMatrixMultiply(X, Y, R);
  R *= -1;
  R.addDiagonal(beta*mu);
}

// R = beta mu I - X Y - dX dY
//
void computeCorrectorRMatrix(const Real &beta,
                             const Real &mu,
                             BlockDiagonalMatrix &X,
                             BlockDiagonalMatrix &dX,
                             BlockDiagonalMatrix &Y,
                             BlockDiagonalMatrix &dY,
                             BlockDiagonalMatrix &R) {
  blockDiagonalMatrixScaleMultiplyAdd(-1, X,  Y,  0, R);
  blockDiagonalMatrixScaleMultiplyAdd(-1, dX, dY, 1, R);
  R.addDiagonal(beta*mu);
}

// Eigenvalues of A, via the QR method
// Inputs:
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
//
void eigenvaluesViaQR(Matrix &A, Vector &workspace, Vector &eigenvalues) {
  assert(A.rows == A.cols);
  assert((int)eigenvalues.size() == A.rows);
  assert((int)workspace.size() == 3*A.rows - 1);

  mpackint info;
  mpackint workSize = workspace.size();
  Rsyev("NoEigenvectors", "LowerTriangular", A.rows, &A.elements[0], A.rows, &eigenvalues[0], &workspace[0], workSize, &info);
  assert(info == 0);
}

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
//
Real minEigenvalueViaQR(Matrix &A, Vector &workspace, Vector &eigenvalues) {
  eigenvaluesViaQR(A, workspace, eigenvalues);
  return eigenvalues[0];
}

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric BlockDiagonalMatrix
// eigenvalues : vector<Vector> of length A.blocks.size()
// workspace   : vector<Vector> of length A.blocks.size()
//
Real minEigenvalueViaQR(BlockDiagonalMatrix &A, vector<Vector> &workspace, vector<Vector> &eigenvalues) {
  assert(A.blocks.size() == eigenvalues.size());
  assert(A.blocks.size() == workspace.size());

  // TODO: get rid of this hack
  Real lambdaMin = 1e100;
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++) {
    Real minBlockLambda = minEigenvalueViaQR(A.blocks[b], workspace[b], eigenvalues[b]);
    #pragma omp critical
    {
      lambdaMin = min(lambdaMin, minBlockLambda);
    }
  }

  return lambdaMin;
}

Real stepLength(BlockDiagonalMatrix &XCholesky,
                BlockDiagonalMatrix &dX,
                BlockDiagonalMatrix &XInvDX,
                vector<Vector> &eigenvalues,
                vector<Vector> &workspace,
                const SDPSolverParameters &parameters) {

  // XInvDX = L^{-1} dX L^{-1 T}, where X = L L^T
  XInvDX.copyFrom(dX);
  lowerTriangularInverseCongruence(XInvDX, XCholesky);

  const Real lambda = minEigenvalueViaQR(XInvDX, workspace, eigenvalues);
  const Real gamma  = parameters.stepLengthReduction;
  if (lambda > -gamma)
    return 1;
  else
    return -gamma/lambda;
}

void computeSchurComplementCholesky(const SDP &sdp,
                                    const BlockDiagonalMatrix &BilinearPairingsXInv,
                                    const BlockDiagonalMatrix &BilinearPairingsY,
                                    const Matrix &BasicKernelSpan,
                                    BlockDiagonalMatrix &SchurBlocks,
                                    BlockDiagonalMatrix &SchurBlocksCholesky,
                                    Matrix &SchurUpdateLowRank,
                                    Matrix &Q,
                                    Matrix &M) {

  timers.computeSchurBlocks.resume();
  computeSchurBlocks(sdp, BilinearPairingsXInv, BilinearPairingsY, SchurBlocks);
  timers.computeSchurBlocks.stop();

  timers.schurBlocksCholesky.resume();
  choleskyDecomposition(SchurBlocks, SchurBlocksCholesky);
  timers.schurBlocksCholesky.stop();

  // SchurUpdateLowRank = (-1 \\ E)
  SchurUpdateLowRank.copyFrom(BasicKernelSpan);

  // SchurUpdateLowRank = SchurBlocksCholesky^{-1} ( - 1 \\ E)
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, SchurUpdateLowRank);

  // Q = SchurUpdateLowRank^T SchurUpdateLowRank
  matrixSquare(SchurUpdateLowRank, Q);

  // Q = M M^T
  choleskyDecomposition(Q, M);

  // SchurUpdateLowRank = SchurUpdateLowRank M^{-T}
  Rtrsm("Right", "Lower", "Transpose", "NonUnitDiagonal",
        SchurUpdateLowRank.rows, SchurUpdateLowRank.cols, 1,
        &M.elements[0], M.rows,
        &SchurUpdateLowRank.elements[0], SchurUpdateLowRank.rows);
}

void printInfoHeader() {
  cout << "     mu       P-obj       D-obj     gap         P-err        D-err       P-step   D-step   beta\n";
  cout << "---------------------------------------------------------------------------------------------------\n";
}

void printInfo(int iteration,
               Real mu,
               SDPSolverStatus status,
               bool isPrimalFeasible,
               bool isDualFeasible,
               Real primalStepLength,
               Real dualStepLength,
               Real betaCorrector) {
  gmp_fprintf(stdout,
              "%3d  %4.1Fe  %+7.2Fe  %+7.2Fe  %+7.2Fe  %s%+7.2Fe%s  %s%+7.2Fe%s  %4.1Fe  %4.1Fe  %4.2Fe\n",
              iteration,
              mu.get_mpf_t(),
              status.primalObjective.get_mpf_t(),
              status.dualObjective.get_mpf_t(),
              status.dualityGap().get_mpf_t(),
              isPrimalFeasible ? "|" : " ", status.primalError.get_mpf_t(), isPrimalFeasible ? "|" : " ",
              isDualFeasible   ? "|" : " ", status.dualError.get_mpf_t(),   isDualFeasible   ? "|" : " ",
              primalStepLength.get_mpf_t(),
              dualStepLength.get_mpf_t(),
              betaCorrector.get_mpf_t());
}

void solveSchurComplementEquation(BlockDiagonalMatrix &SchurBlocksCholesky,
                                  Matrix &SchurUpdateLowRank,
                                  Vector &k,
                                  Vector &dx) {
  // dx = SchurBlocksCholesky^{-1} dx
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, dx);

  // k = -SchurUpdateLowRank^T dx
  vectorScaleMatrixMultiplyTransposeAdd(-1, SchurUpdateLowRank, dx, 0, k);

  // dx = dx + SchurUpdateLowRank k
  vectorScaleMatrixMultiplyAdd(1, SchurUpdateLowRank, k, 1, dx);

  // dx = SchurBlocksCholesky^{-T} dx
  blockMatrixLowerTriangularTransposeSolve(SchurBlocksCholesky, dx);
}

void SDPSolver::computeSearchDirectionWithRMatrix(const BlockDiagonalMatrix &R,
                                                  const bool isPrimalFeasible) {

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
  blockDiagonalMatrixMultiply(PrimalResidues, Y, Z);
  Z -= R;
  blockMatrixSolveWithCholesky(XCholesky, Z);
  Z.symmetrize();

  // dx_k = -d_k + Tr(F_k Z)
  computeSchurRHS(sdp, dualResidues, Z, dx);

  // dx_N = B_{NN}^{-1} dx_N
  // dx_B = E^T dx_N
  solveSchurComplementEquation(SchurBlocksCholesky,
                               SchurUpdateLowRank,
                               basicKernelCoords, dx);

  // dX = R_p + sum_p F_p dx_p
  constraintMatrixWeightedSum(sdp, dx, dX);
  dX += PrimalResidues;
  
  // dY = Symmetrize(X^{-1} (R - dX Y))
  blockDiagonalMatrixMultiply(dX, Y, dY);
  dY -= R;
  blockMatrixSolveWithCholesky(XCholesky, dY);
  dY.symmetrize();
  dY *= -1;
}

SDPSolverStatus SDPSolver::run(const SDPSolverParameters &parameters,
                               const path outFile,
                               const path checkpointFile) {
  printInfoHeader();
  timers.runSolver.resume();

  for (int iteration = 1; iteration <= parameters.maxIterations; iteration++) {
    // Maintain the invariant x_B = g + E^T x_N
    basicCompletion(g, E, x);

    choleskyDecomposition(X, XCholesky);
    choleskyDecomposition(Y, YCholesky);

    timers.bilinearPairings.resume();
    computeInvBilinearPairingsWithCholesky(XCholesky, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsXInv);
    computeBilinearPairings(Y, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsY);
    timers.bilinearPairings.stop();

    // d_k = c_k - Tr(F_k Y)
    computeDualResidues(sdp, Y, BilinearPairingsY, dualResidues);
    nonBasicShift(E, dualResidues, dualResiduesTilde);

    // PrimalResidues = sum_p F_p x_p - X - F_0 (F_0 is zero for now)
    computePrimalResidues(sdp, x, X, PrimalResidues);

    status.primalError     = PrimalResidues.maxAbsElement();
    status.dualError       = maxAbsVectorElement(dualResiduesTilde);
    status.primalObjective = primalObjective(sdp, x);
    status.dualObjective   = dualObjective(sdp, g, dualResidues);

    const bool isPrimalFeasible = status.isPrimalFeasible(parameters);
    const bool isDualFeasible   = status.isDualFeasible(parameters);
    const bool isOptimal        = status.isOptimal(parameters);
    const bool reductionSwitch  = true;

    if (isPrimalFeasible && isDualFeasible && isOptimal) break;

    timers.schurComplementCholesky.resume();
    computeSchurComplementCholesky(sdp,
                                   BilinearPairingsXInv,
                                   BilinearPairingsY,
                                   BasicKernelSpan,
                                   SchurBlocks,
                                   SchurBlocksCholesky,
                                   SchurUpdateLowRank,
                                   Q,
                                   M);
    timers.schurComplementCholesky.stop();

    Real mu = frobeniusProductSymmetric(X, Y)/X.dim;

    // Mehrotra predictor solution for (dx, dX, dY)
    Real betaPredictor = predictorCenteringParameter(parameters, reductionSwitch,
                                                     isPrimalFeasible && isDualFeasible);
    timers.predictorSolution.resume();
    computePredictorRMatrix(betaPredictor, mu, X, Y, R);
    computeSearchDirectionWithRMatrix(R, isPrimalFeasible);
    timers.predictorSolution.stop();

    // Mehrotra corrector solution for (dx, dX, dY)
    Real betaCorrector = correctorCenteringParameter(parameters, X, dX, Y, dY, mu,
                                                     isPrimalFeasible && isDualFeasible);
    timers.correctorSolution.resume();
    computeCorrectorRMatrix(betaCorrector, mu, X, dX, Y, dY, R);
    computeSearchDirectionWithRMatrix(R, isPrimalFeasible);
    timers.correctorSolution.stop();

    // Step length to preserve positive definiteness
    Real primalStepLength = stepLength(XCholesky, dX, StepMatrixWorkspace,
                                       eigenvaluesWorkspace, QRWorkspace, parameters);
    Real dualStepLength   = stepLength(YCholesky, dY, StepMatrixWorkspace,
                                       eigenvaluesWorkspace, QRWorkspace, parameters);

    // Update current point
    vectorScaleMultiplyAdd(primalStepLength, dx, 1, x);
    dX *= primalStepLength;
    X += dX;
    dY *= dualStepLength;
    Y += dY;

    printInfo(iteration, mu, status, isPrimalFeasible, isDualFeasible,
              primalStepLength, dualStepLength, betaCorrector);
  }

  timers.runSolver.stop();
  return status;
}

void printSDPDenseFormat(ostream& os, const SDP &sdp) {
  BlockDiagonalMatrix F(BlockDiagonalMatrix(sdp.psdMatrixBlockDims()));

  os << "* SDP dense format" << endl;
  os << sdp.affineConstants.size() << " = mDIM" << endl;
  os << F.blocks.size() + 1 << " = nBLOCK" << endl;
  os << "{";
  for (unsigned int b = 0; b < F.blocks.size(); b++) {
    os << F.blocks[b].rows;
    if (b != F.blocks.size() - 1)
      os << ", ";
  }
  os << "} = bLOCKsTRUCT" << endl;

  os << sdp.affineConstants << endl;

  F *= 0;
  os << F << endl;

  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t++) {
      F *= 0;

      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end();
           b++) {
        const int delta = sdp.bilinearBases[*b].rows;
        const Real *q = &sdp.bilinearBases[*b].elements[(t->k) * delta];

        for (int e = 0; e < delta; e++) {
          for (int f = 0; f < delta; f++) {
            F.blocks[*b].elt((t->r)*delta + e, (t->s)*delta + f) = (*(q+e)) * (*(q+f));
          }
        }
        F.blocks[*b].symmetrize();
      }

      os << F << endl;
    }
  }

}

void printSDPBHeader(const path &sdpFile,
                     const path &outFile,
                     const path &checkpointFile,
                     const SDPSolverParameters &parameters) {
  cout << "SDPB started at " << second_clock::local_time() << endl;
  cout << "SDP file        : " << sdpFile        << endl;
  cout << "out file        : " << outFile        << endl;
  cout << "checkpoint file : " << checkpointFile << endl;
  cout << "using " << omp_get_max_threads() << " threads." << endl;

  cout << "\nParameters:\n";
  cout << parameters << endl;
}

void solveSDP(const path sdpFile,
              const optional<path> outFileOption,
              const optional<path> checkpointFileOption,
              const optional<path> paramFileOption) {

  path outFile = sdpFile;
  outFile.replace_extension("out");
  if (outFileOption)
    outFile = *outFileOption;

  path checkpointFile = sdpFile;
  checkpointFile.replace_extension("ck");
  if (checkpointFileOption)
    checkpointFile = *checkpointFileOption;
    
  SDPSolverParameters parameters;
  if (paramFileOption)
    parameters = SDPSolverParameters(*paramFileOption);

  mpf_set_default_prec(parameters.precision);
  cout.precision(int(parameters.precision * 0.30102999566398114 + 5));
  omp_set_num_threads(parameters.maxThreads);

  printSDPBHeader(sdpFile, outFile, checkpointFile, parameters);

  const SDP sdp = readBootstrapSDP(sdpFile);

  SDPSolver solver(sdp);
  solver.initialize(parameters);
  SDPSolverStatus status = solver.run(parameters, outFile, checkpointFile);
  
  cout << "\nStatus:\n";
  cout << status << endl;
  cout << timers << endl;

  cout << "X = " << solver.X << ";\n";
  cout << "Y = " << solver.Y << ";\n";
  cout << "x = " << solver.x << ";\n";
  // cout << "BilinearPairingsXInv = " << solver.BilinearPairingsXInv << endl;
  // cout << "BilinearPairingsY = " << solver.BilinearPairingsY << endl;
  // cout << "schurComplement = " << solver.schurComplement << ";\n";
  // cout << "R = " << solver.R << ";\n";
  // cout << "dualResidues = " << solver.dualResidues << ";\n";
  // cout << "PrimalResidues = " << solver.PrimalResidues << ";\n";
  // cout << "Z = " << solver.Z << ";\n";
  // cout << "dx = " << solver.dx << ";\n";
  // cout << "dX = " << solver.dX << ";\n";
  // cout << "dY = " << solver.dY << ";\n";

  // path datFile = sdpFile;
  // datFile.replace_extension("dat");
  // ofstream datStream;
  // datStream.open(datFile.c_str());
  // datStream.precision(parameters.precision);
  //cout << sdp;
  //printSDPDenseFormat(cout, sdp);
  // datStream.close();
}

void testBilinearPairings(const path sdpFile) {
  SDPSolverParameters parameters;
  const SDP sdp = readBootstrapSDP(sdpFile);
  SDPSolver solver(sdp);
  solver.initialize(parameters);
  for (int i = 0; i < 100; i++) {
    cout << "i = " << i << endl;
    computeBilinearPairings(solver.Y, sdp.bilinearBases, solver.bilinearPairingsWorkspace, solver.BilinearPairingsY);
  }
}

void testCholeskyUpdate() {
  Matrix A(4,4);
  Matrix B(A);
  Matrix C(A);
  Matrix L(A);
  Matrix LT(L);
  Matrix V(4, 2);
  Matrix VT(V.cols, V.rows);
  V.elt(0,0) =1;
  V.elt(1,0) =2;
  V.elt(2,0) =3;
  V.elt(3,0) =4;
  V.elt(0,1) =5;
  V.elt(1,1) =4;
  V.elt(2,1) =3;
  V.elt(3,1) =2;
  for (int r = 0; r < V.rows; r++)
    for (int c = 0; c < V.cols; c++)
      VT.elt(c, r) = V.elt(r,c);
  Matrix U(V);

  A.addDiagonal(4);
  cout << "A = " << A << endl;
  cout << "V = " << V << endl;
  choleskyDecomposition(A, L);
  choleskyUpdate(L, U);
  transpose(L, LT);

  matrixMultiply(V, VT, B);
  B += A;
  matrixMultiply(L, LT, C);
  C -= B;

  cout << "L L^T - (A + V V^T) = " << C << endl;
}

void testMatrix() {
  Matrix A(3,3);
  A.elt(0,0) = 1;
  A.elt(1,0) = 2;
  A.elt(2,0) = 3;
  A.symmetrize();
  cout << A << endl;
}

const char *help(const char *cmd) {
  return cmd;
}

int main(int argc, char** argv) {

  path sdpFile;
  optional<path> outFile;
  optional<path> checkpointFile;
  optional<path> paramFile;

  int i = 1;
  while (i < argc) {
    const char *arg = argv[i];
    if (!strcmp(arg, "-h") || !strcmp(arg, "--help")) {
      cout << help(argv[0]);
      return 0;
    } else if ((!strcmp(arg, "-s") || !strcmp(arg, "--sdpFile"))        && i+1 < argc) {
      sdpFile = path(argv[i+1]);
      i += 2;
    } else if ((!strcmp(arg, "-p") || !strcmp(arg, "--paramFile"))      && i+1 < argc) {
      paramFile = path(argv[i+1]);
      i += 2;
    } else if ((!strcmp(arg, "-o") || !strcmp(arg, "--outFile"))        && i+1 < argc) {
      outFile = path(argv[i+1]);
      i += 2;
    } else if ((!strcmp(arg, "-c") || !strcmp(arg, "--checkpointFile")) && i+1 < argc) {
      checkpointFile = path(argv[i+1]);
      i += 2;
    } else {
      cout << help(argv[0]);
      return 1;
    }
  }

  solveSDP(sdpFile, outFile, checkpointFile, paramFile);
  //testMatrix();
  //testBilinearPairings(sdpFile);

  //testBlockCongruence();
  //testBlockDiagonalCholesky();
  //testSDPSolver(argv[1], argv[2]);
  //testCholeskyUpdate();
  //testMinEigenvalue();
  //testTensorCongruence();
  return 0;
}
