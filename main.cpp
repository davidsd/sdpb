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
  if (last > 0)
    os << v[last];
  os << "}";
  return os;
}

typedef vector<Real> Vector;

Real vectorTotal(const Vector &v) {
  Real total = 0;
  for (Vector::const_iterator x = v.begin(); x != v.end(); x++)
    total += *x;
  return total;
}

Real maxAbsVectorElement(const Vector &v) {
  Real max = abs(v[0]);
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

  void transposeInplace() {
    assert (rows == cols);
    for (int c = 0; c < cols; c++) {
      for (int r = 0; r < c; r++) {
        Real tmp = elt(r, c);
        elt(r, c) = elt(c, r);
        elt(c, r) = tmp;
      }
    }
  }

  void addColumn() {
    cols += 1;
    elements.resize(rows*cols);
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

Real dotProduct(const Vector &u, const Vector v) {
  Real result = 0;
  for (unsigned int i = 0; i < u.size(); i++)
    result += u[i]*v[i];
  return result;
}

Real norm(const Vector &v) {
  return sqrt(dotProduct(v,v));
}

// y := alpha*A*x + beta*y
//
void vectorScaleMatrixMultiplyAdd(Real alpha, Matrix &A, Vector &x, Real beta, Vector &y) {
  assert(A.cols == (int)x.size());
  assert(A.rows == (int)y.size());

  Rgemv("NoTranspose",
        A.rows, A.cols, alpha,
        &A.elements[0], A.rows,
        &x[0], (int)x.size(),
        beta,
        &y[0], (int)y.size());
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
  cout << info << endl;
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

// result = A^-1
// Inputs:
// - A      : dim x dim lower-triangular matrix
// - result : dim x dim lower-triangular matrix
//
void inverseLowerTriangular(Matrix &L, Matrix &result) {
  int dim = L.rows;
  assert(result.rows == dim);
  assert(result.cols == dim);

  result.setIdentity();
  lowerTriangularSolve(L, &result.elements[0], dim, dim);
}

// result = choleskyDecomposition(A)^-1
// Inputs:
// - A      : dim x dim symmetric matrix
// - work   : dim x dim matrix
// - result : dim x dim lower-triangular matrix
//
void inverseCholesky(Matrix &A, Matrix &work, Matrix &result) {
  choleskyDecomposition(A, work);
  inverseLowerTriangular(work, result);
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

// invCholesky = choleskyDecomposition(a)^-1
// inverse = a^-1
// Inputs:
// - a           : dim x dim symmetric matrix
// - work        : dim x dim matrix
// - invCholesky : dim x dim lower-triangular matrix
// - inverse     : dim x dim symmetric matrix
//
// TODO: we can save memory by using inverse as the work matrix for
// inverse cholesky
//
void inverseCholeskyAndInverse(Matrix &a, Matrix &work, Matrix &invCholesky, Matrix &inverse) {
  int dim = a.rows;
  assert(a.cols == dim);
  assert(work.rows        == dim && work.cols        == dim);
  assert(invCholesky.rows == dim && invCholesky.cols == dim);
  assert(inverse.rows     == dim && inverse.cols     == dim);

  inverseCholesky(a, work, invCholesky);

  inverse.elements = invCholesky.elements;
  Rtrmm("Left", "Lower", "Transpose", "NonUnitDiag", dim, dim, 1,
        &invCholesky.elements[0], dim,
        &inverse.elements[0], dim);
}

// X := AInvCholesky^T AInvCholesky X
// Inputs:
// - AInvCholesky : dim x dim lower triangular matrix
// - X            : dim x dim matrix
//
void matrixSolveWithInverseCholesky(Matrix &AInvCholesky, Matrix &X) {
  int dim = X.rows;
  assert(X.cols == dim);
  assert(AInvCholesky.rows == dim);
  assert(AInvCholesky.cols == dim);

  Rtrmm("Left", "Lower", "NoTranspose", "NonUnitDiag", dim, dim, 1,
        &AInvCholesky.elements[0], dim,
        &X.elements[0], dim);

  Rtrmm("Left", "Lower", "Transpose", "NonUnitDiag", dim, dim, 1,
        &AInvCholesky.elements[0], dim,
        &X.elements[0], dim);
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



void testTensorCongruence() {
  Matrix a(4,4);
  Matrix b(2,3);
  Matrix result(6,6);
  Matrix work(4,6);
  a.setIdentity();
  b.elt(0,0) =2;
  b.elt(1,0) =3;
  b.elt(0,1) =4;
  b.elt(1,1) =5;
  b.elt(0,2) =6;
  b.elt(1,2) =7;

  tensorMatrixCongruenceTranspose(a, b, work, result);

  cout << a << endl;
  cout << b << endl;
  cout << work << endl;
  cout << result << endl;
  
}

class BlockDiagonalMatrix {
public:
  int dim;
  Vector diagonalPart;
  vector<Matrix> blocks;
  vector<int> blockStartIndices;

  BlockDiagonalMatrix(int diagonalSize, const vector<int> &blockSizes):
    dim(diagonalSize),
    diagonalPart(Vector(diagonalSize, 0)) {

    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(Matrix(blockSizes[i], blockSizes[i]));
      blockStartIndices.push_back(dim);
      dim += blockSizes[i];
    }
  }

  void setZero() {
    fillVector(diagonalPart, 0);

    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].setZero();
  }

  void addDiagonal(const Real &c) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] += c;

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].addDiagonal(c);
  }

  void setIdentity() {
    setZero();
    addDiagonal(1);
  }

  void addDiagonalPart(const Vector &v, const Real &alpha) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] += alpha*v[i];
  }

  void operator+=(const BlockDiagonalMatrix &A) {
    addDiagonalPart(A.diagonalPart, 1);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] += A.blocks[b];
  }

  void operator-=(const BlockDiagonalMatrix &A) {
    addDiagonalPart(A.diagonalPart, -1);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] -= A.blocks[b];
  }

  void operator*=(const Real &c) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] *= c;

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] *= c;
  }

  void copyFrom(const BlockDiagonalMatrix &A) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] = A.diagonalPart[i];
    
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].copyFrom(A.blocks[b]);
  }

  // fill out a BlockDiagonalMatrix into a full Matrix A
  void copyInto(Matrix &A) {
    A.setZero();

    int p = 0;
    for(; p < (int)diagonalPart.size(); p++)
      A.elt(p, p) = diagonalPart[p];

    // TODO: parallelize this loop
    for (unsigned int b = 0; b < blocks.size(); b++) {
      int dim = blocks[b].cols;
      for (int c = 0; c < dim; c++)
        for (int r = 0; r < dim; r++)
          A.elt(p + r, p + c) = blocks[b].elt(r, c);
      p += dim;
    }
  }

  void symmetrize() {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].symmetrize();
  }

  Real maxAbsElement() const {
    Real max = maxAbsVectorElement(diagonalPart);
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

  os << "{{";
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++) {
    os << A.diagonalPart[i];
    if (i != A.diagonalPart.size() - 1)
       os << ", ";
  }
  os << "}, {";

  for (unsigned int b = 0; b < A.blocks.size(); b++) {
    os << A.blocks[b];
    if (b < A.blocks.size() - 1)
      os << ", ";
  }
  os << "}}";
  return os;
}

Real frobeniusProductSymmetric(const BlockDiagonalMatrix &A,
                               const BlockDiagonalMatrix &B) {
  Real result = dotProduct(A.diagonalPart, B.diagonalPart);

  #pragma omp parallel for schedule(dynamic)
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
  for (unsigned int i = 0; i < X.diagonalPart.size(); i++)
    result += (X.diagonalPart[i] + dX.diagonalPart[i])*(Y.diagonalPart[i] + dY.diagonalPart[i]);

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < X.blocks.size(); b++)
    result += frobeniusProductOfSums(X.blocks[b], dX.blocks[b], Y.blocks[b], dY.blocks[b]);

  return result;
}

// void blockDiagonalMatrixAdd(const BlockDiagonalMatrix &A,
//                             const BlockDiagonalMatrix &B,
//                             BlockDiagonalMatrix &result) {
//   for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
//     result.diagonalPart[i] = A.diagonalPart[i] + B.diagonalPart[i];

//   for (unsigned int b = 0; b < A.blocks.size(); b++) 
//     matrixAdd(A.blocks[b], B.blocks[b], result.blocks[b]);
// }

void blockDiagonalMatrixScaleMultiplyAdd(Real alpha,
                                         BlockDiagonalMatrix &A,
                                         BlockDiagonalMatrix &B,
                                         Real beta,
                                         BlockDiagonalMatrix &C) {
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
    C.diagonalPart[i] = alpha*A.diagonalPart[i]*B.diagonalPart[i] + beta*C.diagonalPart[i];

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    matrixScaleMultiplyAdd(alpha, A.blocks[b], B.blocks[b], beta, C.blocks[b]);
}

void blockDiagonalMatrixMultiply(BlockDiagonalMatrix &A,
                                 BlockDiagonalMatrix &B,
                                 BlockDiagonalMatrix &C) {
  blockDiagonalMatrixScaleMultiplyAdd(1, A, B, 0, C);
}

void lowerTriangularCongruence(BlockDiagonalMatrix &A, BlockDiagonalMatrix &L) {
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
    A.diagonalPart[i] *= L.diagonalPart[i]*L.diagonalPart[i];

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    lowerTriangularCongruence(A.blocks[b], L.blocks[b]);
}

void choleskyDecomposition(BlockDiagonalMatrix &A,
                           BlockDiagonalMatrix &L) {
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
    L.diagonalPart[i] = sqrt(A.diagonalPart[i]);

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    choleskyDecomposition(A.blocks[b], L.blocks[b]);
}

void inverseCholesky(BlockDiagonalMatrix &A,
                     BlockDiagonalMatrix &work,
                     BlockDiagonalMatrix &AInvCholesky) {
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
    AInvCholesky.diagonalPart[i] = 1/sqrt(A.diagonalPart[i]);

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    inverseCholesky(A.blocks[b], work.blocks[b], AInvCholesky.blocks[b]);
}

void inverseCholeskyAndInverse(BlockDiagonalMatrix &A,
                               BlockDiagonalMatrix &work,
                               BlockDiagonalMatrix &AInvCholesky,
                               BlockDiagonalMatrix &AInv) {

  for (unsigned int i = 0; i < A.diagonalPart.size(); i++) {
    Real d = A.diagonalPart[i];
    AInvCholesky.diagonalPart[i] = 1/sqrt(d);
    AInv.diagonalPart[i] = 1/d;
  }

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < A.blocks.size(); b++)
    inverseCholeskyAndInverse(A.blocks[b], work.blocks[b], AInvCholesky.blocks[b], AInv.blocks[b]);
}

void blockMatrixSolveWithInverseCholesky(BlockDiagonalMatrix &AInvCholesky,
                                         BlockDiagonalMatrix &X) {
  for (unsigned int i = 0; i < X.diagonalPart.size(); i++)
    X.diagonalPart[i] *= AInvCholesky.diagonalPart[i] * AInvCholesky.diagonalPart[i];

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < X.blocks.size(); b++)
    matrixSolveWithInverseCholesky(AInvCholesky.blocks[b], X.blocks[b]);
}

void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, Matrix &B) {
  // TODO: allow diagonal parts for completeness
  assert(L.diagonalPart.size() == 0);

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularSolve(L.blocks[b], &B.elt(L.blockStartIndices[b], 0), B.cols, B.rows);
}

void blockMatrixLowerTriangularSolve(BlockDiagonalMatrix &L, Vector &v) {
  // TODO: allow diagonal parts for completeness
  assert(L.diagonalPart.size() == 0);

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularSolve(L.blocks[b], &v[L.blockStartIndices[b]], 1, v.size());
}

void blockMatrixLowerTriangularTransposeSolve(BlockDiagonalMatrix &L, Vector &v) {
  // TODO: allow diagonal parts for completeness
  assert(L.diagonalPart.size() == 0);

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < L.blocks.size(); b++)
    lowerTriangularTransposeSolve(L.blocks[b], &v[L.blockStartIndices[b]], 1, v.size());
}

void testBlockDiagonalCholesky() {
  vector<int> sizes;
  sizes.push_back(3);
  sizes.push_back(4);

  BlockDiagonalMatrix a(2, sizes);
  a.setIdentity();
  a.diagonalPart[0] = 2;
  a.diagonalPart[1] = 3;
  Real aBlock0[] = {14,3,8,3,10,9,8,9,14};
  a.blocks[0].elements.insert(a.blocks[0].elements.begin(), aBlock0, aBlock0 + 9);

  BlockDiagonalMatrix work(2, sizes);
  BlockDiagonalMatrix invCholesky(2, sizes);
  BlockDiagonalMatrix inverse(2, sizes);

  inverseCholeskyAndInverse(a, work, invCholesky, inverse);

  cout << a << endl;
  cout << invCholesky << endl;
  cout << inverse << endl;
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

  void addSlack() {
    polMatrixValues.addColumn();
    for (int r = 0; r < polMatrixValues.rows; r++) {
      Real total = 0;
      for (int c = 0; c < polMatrixValues.cols - 1; c++)
        total += polMatrixValues.elt(r, c);
      polMatrixValues.elt(r, polMatrixValues.cols - 1) = -total;
    }
    objective.push_back(-vectorTotal(objective));
  }

  void rescale() {
    Vector colWeightedNormSq(objective.size());
    
    for (int p = 0; p < polMatrixValues.rows; p++) {
      Real rowNormSq = 0;
      for (int n = 0; n < polMatrixValues.cols; n++)
        rowNormSq += polMatrixValues.elt(p, n)*polMatrixValues.elt(p, n);

      for (int n = 0; n < polMatrixValues.cols; n++)
        colWeightedNormSq[n] += polMatrixValues.elt(p, n)*polMatrixValues.elt(p, n) / rowNormSq;
    }

    Vector rescaling(objective.size());
    for (unsigned int n = 0; n < rescaling.size(); n++) {
      rescaling[n] = polMatrixValues.rows/sqrt(colWeightedNormSq[n]);
      objective[n] *= rescaling[n];
      for (int p = 0; p < polMatrixValues.rows; p++)
        polMatrixValues.elt(p, n) *= rescaling[n];
    }
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

class Polynomial {
public:
  Vector coeffs;

  Polynomial(): coeffs(1, 0) {}

  int degree() const {
    return coeffs.size() - 1;
  };

  Real operator()(const Real &x) const {
    int deg = degree();
    Real y = coeffs[deg];
    for (int i = deg - 1; i >= 0; i--) {
      y *= x;
      y += coeffs[i];
    }
    return y;
  }

  friend ostream& operator<<(ostream& os, const Polynomial& p);

};

ostream& operator<<(ostream& os, const Polynomial& p) {
  for (int i = p.degree(); i >= 0; i--) {
    os << p.coeffs[i];
    if (i > 1)
      os << "x^" << i << " + ";
    else if (i == 1)
      os << "x + ";
  }
  return os;
}

class PolynomialVectorMatrix {
public:
  int rows;
  int cols;
  vector<vector<Polynomial> > elements;
  
  const vector<Polynomial>& elt(const int r, const int c) const {
    return elements[r+c*rows];
  }

  int degree() const {
    int d = 0;
    for (vector<vector<Polynomial> >::const_iterator e = elements.begin(); e != elements.end(); e++)
      for (vector<Polynomial>::const_iterator p = e->begin(); p != e->end(); p++)
        d = max(p->degree(), d);
    return d;
  }

};

Matrix basisAtSamplePoints(int basisSize, int numPoints, bool withSqrt,
                           const vector<Polynomial> &basisPols,
                           const Vector &samplePoints) {
  Matrix m(basisSize, numPoints);
  for (int k = 0; k < numPoints; k++) {
    Real x = samplePoints[k];
    Real sqrtX = sqrt(x);

    for (int n = 0; n < basisSize; n++) {
      if (withSqrt)
        m.elt(n, k) = basisPols[n](x)*sqrtX;
      else
        m.elt(n, k) = basisPols[n](x);
    }
  }

  return m;
}

SDP bootstrapSDP(const Vector &objective,
                 const Vector &normalization,
                 const vector<PolynomialVectorMatrix> &positiveMatrixPols,
                 const vector<Polynomial> &bilinearBasisPols,
                 const Vector &polynomialSamplePoints) {
  SDP sdp;
  sdp.objective = objective;

  int numConstraints = 0;
  for (vector<PolynomialVectorMatrix>::const_iterator m = positiveMatrixPols.begin();
       m != positiveMatrixPols.end();
       m++) {
    int dimension = m->cols;
    int degree    = m->degree();

    sdp.dimensions.push_back(dimension);
    sdp.degrees.push_back(degree);
    numConstraints += (degree+1)*dimension*(dimension+1)/2;
  }

  // For the normalization constraint
  // sdp.dimensions.push_back(1);
  // sdp.degrees.push_back(0);
  // numConstraints += 1;

  sdp.polMatrixValues = Matrix(numConstraints, sdp.objective.size());
  sdp.affineConstants = Vector(numConstraints, 0);

  // normalization constraint
  // sdp.affineConstants[numConstraints-1] = 1;

  int p = 0;
  for (vector<PolynomialVectorMatrix>::const_iterator m = positiveMatrixPols.begin();
       m != positiveMatrixPols.end();
       m++) {

    int degree = m->degree();
    int delta1 = degree/2;
    int delta2 = (degree-1)/2;

    vector<int> blocks;

    blocks.push_back(sdp.bilinearBases.size());
    sdp.bilinearBases.push_back(basisAtSamplePoints(delta1+1, degree+1, false,
                                                    bilinearBasisPols,
                                                    polynomialSamplePoints));

    if (delta2 >= 0) {
      blocks.push_back(sdp.bilinearBases.size());
      sdp.bilinearBases.push_back(basisAtSamplePoints(delta2+1, degree+1, true,
                                                      bilinearBasisPols,
                                                      polynomialSamplePoints));
    }

    sdp.blocks.push_back(blocks);

    for (int s = 0; s < m->cols; s++) {
      for (int r = 0; r <= s; r++) {
        for (int k = 0; k <= degree; k++, p++) {
          const Real xk = polynomialSamplePoints[k];
          for (unsigned int n = 0; n < sdp.objective.size(); n++)
            sdp.polMatrixValues.elt(p, n) = -m->elt(r,s)[n](xk);
        }
      }
    }
  }
  // assert(p == numConstraints-1);
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

SDP bootstrapSDP2(const Vector &objective,
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
  sdp.affineConstants = Vector(numConstraints, 0);

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


SDP parseBootstrapSDP2(XMLElement *sdpXml) {
  return bootstrapSDP2(parseVector(sdpXml->FirstChildElement("objective")),
                       parseVector(sdpXml->FirstChildElement("normalization")),
                       parseMany("sampledMatrixPolynomial",
                                 parseSampledMatrixPolynomial,
                                 sdpXml->FirstChildElement("sampledPositiveMatrices")));
}

Polynomial parsePolynomial(XMLElement *polXml) {
  Polynomial p;
  p.coeffs = parseMany("coeff", parseReal, polXml);
  return p;
}

vector<Polynomial> parsePolynomialVector(XMLElement *polVecXml) {
  return parseMany("polynomial", parsePolynomial, polVecXml);
}

PolynomialVectorMatrix parsePolynomialVectorMatrix(XMLElement *polVecMatrixXml) {
  PolynomialVectorMatrix m;
  m.rows = parseInt(polVecMatrixXml->FirstChildElement("rows"));
  m.cols = parseInt(polVecMatrixXml->FirstChildElement("cols"));
  m.elements = parseMany("polynomialVector", parsePolynomialVector,
                         polVecMatrixXml->FirstChildElement("elements"));
  return m;
}

SDP parseBootstrapSDP(XMLElement *sdpXml) {
  return bootstrapSDP(parseVector(sdpXml->FirstChildElement("objective")->FirstChildElement("vector")),
                      parseVector(sdpXml->FirstChildElement("normalization")->FirstChildElement("vector")),
                      parseMany("polynomialVectorMatrix",
                                parsePolynomialVectorMatrix,
                                sdpXml->FirstChildElement("positiveMatrixPols")),
                      parsePolynomialVector(sdpXml->FirstChildElement("bilinearBasisPols")->FirstChildElement("polynomialVector")),
                      parseVector(sdpXml->FirstChildElement("polynomialSamplePoints")->FirstChildElement("vector")));
}

SDP readBootstrapSDP(const path sdpFile) {
  XMLDocument doc;
  doc.LoadFile(sdpFile.c_str());
  return parseBootstrapSDP2(doc.FirstChildElement("sdp"));
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
    maxIterations(500),
    dualityGapThreshold("1e-100"),
    primalFeasibilityThreshold("1e-30"),
    dualFeasibilityThreshold("1e-30"),
    initialMatrixScale("1e20"),
    feasibleCenteringParameter("0.1"),
    infeasibleCenteringParameter("0.3"),
    stepLengthReduction("0.7"),
    precision(400),
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
  BlockDiagonalMatrix PrimalResidues;

  // Schur complement for computing search direction
  Matrix SchurComplementCholesky;

  // For free variable elimination
  Matrix E;
  Matrix ET;
  Vector g;
  Vector y;
  Vector dualResiduesTilde;
  Real c0Tilde;
  //  Vector cTilde;

  // New variables for Schur complement calculation
  // Matrix schurComplementP;
  // Matrix schurComplementQ;
  // Vector schurComplementY;
  // Vector schurComplementWork;

  // intermediate computations
  BlockDiagonalMatrix XInv;
  BlockDiagonalMatrix XInvCholesky;
  BlockDiagonalMatrix YInvCholesky;
  BlockDiagonalMatrix Z;
  BlockDiagonalMatrix R;
  BlockDiagonalMatrix BilinearPairingsXInv;
  BlockDiagonalMatrix BilinearPairingsY;
  BlockDiagonalMatrix SchurBlocks;
  BlockDiagonalMatrix SchurBlocksCholesky;
  Matrix SchurUpdateLowRank;

  // additional workspace variables
  BlockDiagonalMatrix XInvWorkspace;
  BlockDiagonalMatrix StepMatrixWorkspace;
  vector<Matrix> bilinearPairingsWorkspace;
  vector<Vector> eigenvaluesWorkspace;
  vector<Vector> QRWorkspace;

  SDPSolver(const SDP &sdp):
    sdp(sdp),
    x(sdp.numConstraints(), 0),
    X(0, sdp.psdMatrixBlockDims()),
    Y(X),
    dx(x),
    dX(X),
    dY(Y),
    dualResidues(x),
    PrimalResidues(X),
    SchurComplementCholesky(sdp.numConstraints(), sdp.numConstraints()),
    E(sdp.numConstraints() - sdp.objective.size(), sdp.objective.size()),
    ET(E.cols, E.rows),
    g(sdp.objective.size()),
    y(x),
    dualResiduesTilde(sdp.numConstraints() - sdp.objective.size()),
    // cTilde(sdp.numConstraints() - sdp.objective.size()),
    // schurComplementP(sdp.polMatrixValues.cols, sdp.polMatrixValues.cols),
    // schurComplementQ(schurComplementP),
    // schurComplementY(sdp.polMatrixValues.rows),
    // schurComplementWork(schurComplementQ.rows),
    XInv(X),
    XInvCholesky(X),
    YInvCholesky(X),
    Z(X),
    R(X),
    BilinearPairingsXInv(0, sdp.bilinearPairingBlockDims()),
    BilinearPairingsY(BilinearPairingsXInv),
    SchurBlocks(0, sdp.schurBlockDims()),
    SchurBlocksCholesky(SchurBlocks),
    SchurUpdateLowRank(sdp.polMatrixValues),
    XInvWorkspace(X),
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
    for (int p = 0; p < ET.cols; p++)
      for (int n = 0; n < ET.rows; n++)
        ET.elt(n, p) = -sdp.polMatrixValues.elt(p + nonBasicStart, n);
    // ET = D_B^{-1 T} ET = -D_B^{-1 T} D_N^T
    solveWithLUDecompositionTranspose(DBLU, DBLUipiv, &ET.elements[0], ET.cols, ET.rows);
    // E = ET^T
    transpose(ET, E);

    // g = -D_B^{-T} f
    for (unsigned int n = 0; n < g.size(); n++)
      g[n] = -sdp.objective[n];
    solveWithLUDecompositionTranspose(DBLU, DBLUipiv, &g[0], 1, g.size());

    // c0Tilde = - c_B^T g
    c0Tilde = 0;
    for (unsigned int n = 0; n < g.size(); n++)
      c0Tilde -= g[n]*sdp.affineConstants[n];

    // We don't actually need this, since it'll be included implicitly
    // in later calculations
    // cTilde = c_N + E c_B
    // for (unsigned int p = 0; p < cTilde.size(); p++) {
    //   cTilde[p] = sdp.affineConstants[p + nonBasicStart];
    //   for (int n = 0; n < E.cols; n++)
    //     cTilde[p] += E.elt(p, n) * sdp.affineConstants[n];
    // }

    cout << "polMatrixValues = " << sdp.polMatrixValues << ";\n";
    cout << "E = " << E << ";\n";
    // cout << "cTilde = " << cTilde << ";\n";
    cout << "c0Tilde = " << c0Tilde << ";\n";
    
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
           
// result[i] = u[i] * v[i]
//                
void componentProduct(const Vector &u, const Vector &v, Vector &result) {
  for (unsigned int i = 0; i < u.size(); i++)
    result[i] = u[i] * v[i];
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

// xTilde_N = x_N + E x_B
void nonBasicShift(const Matrix &E, const Vector &x, Vector &xTilde) {
  int nonBasicStart = x.size() - xTilde.size();
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

      // Y no longer has a diagonal part
      // for (int n = 0; n < sdp.polMatrixValues.cols; n++)
      //   dualResidues[p] -= Y.diagonalPart[n] * sdp.polMatrixValues.elt(p, n);

      dualResidues[p] += sdp.affineConstants[p];
    }
  }
}

void constraintMatrixWeightedSum(const SDP &sdp, const Vector x, BlockDiagonalMatrix &result)  {

  // for (unsigned int n = 0; n < result.diagonalPart.size(); n++) {
  //   result.diagonalPart[n] = 0;
  //   for (unsigned int p = 0; p < x.size(); p++)
  //     result.diagonalPart[n] += x[p]*sdp.polMatrixValues.elt(p, n);
  // }

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

  for (unsigned int p = 0; p < r.size(); p++) {
    r[p] = -dualResidues[p];
    // for (unsigned int n = 0; n < Z.diagonalPart.size(); n++)
    //   r[p] -= sdp.polMatrixValues.elt(p, n)*Z.diagonalPart[n];
  }

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

// Real dualObjective(const SDP &sdp, const BlockDiagonalMatrix &Y) {
//   return dotProduct(sdp.objective, Y.diagonalPart);
// }

Real dualObjective(const SDP &sdp, const Vector &g, const Vector &dualResidues) {
  Real result = 0;
  for (unsigned int i = 0; i < g.size(); i++)
    result += g[i]*(sdp.affineConstants[i] - dualResidues[i]);
  return result;
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
  // Real lambdaMin = A.diagonalPart[0];
  // for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
  //   lambdaMin = min(lambdaMin, A.diagonalPart[i]);

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

Real stepLength(BlockDiagonalMatrix &XInvCholesky,
                BlockDiagonalMatrix &dX,
                BlockDiagonalMatrix &XInvDX,
                vector<Vector> &eigenvalues,
                vector<Vector> &workspace,
                const SDPSolverParameters &parameters) {

  // XInvDX = L^{-1} dX L^{-1 T}, where X = L L^T
  XInvDX.copyFrom(dX);
  lowerTriangularCongruence(XInvDX, XInvCholesky);

  const Real lambda = minEigenvalueViaQR(XInvDX, workspace, eigenvalues);
  const Real gamma  = parameters.stepLengthReduction;
  if (lambda > -gamma)
    return 1;
  else
    return -gamma/lambda;
}

void computeSchurComplementCholesky(const SDP &sdp,
                                    const BlockDiagonalMatrix &XInv,
                                    const BlockDiagonalMatrix &BilinearPairingsXInv,
                                    const BlockDiagonalMatrix &Y,
                                    const BlockDiagonalMatrix &BilinearPairingsY,
                                    BlockDiagonalMatrix &SchurBlocks,
                                    BlockDiagonalMatrix &SchurBlocksCholesky,
                                    Matrix &SchurUpdateLowRank,
                                    Matrix &SchurComplementCholesky,
                                    int i) {

  timers.computeSchurBlocks.resume();
  computeSchurBlocks(sdp, BilinearPairingsXInv, BilinearPairingsY, SchurBlocks);
  timers.computeSchurBlocks.stop();

  // cout << "XInvDiagonal[" << i << "] = " << XInv.diagonalPart << ";\n";
  // cout << "YDiagonal[" << i << "] = " << Y.diagonalPart << ";\n";
  //  cout << "SchurBlocks[" << i << "] = " << SchurBlocks << ";\n";

  timers.schurBlocksCholesky.resume();
  choleskyDecomposition(SchurBlocks, SchurBlocksCholesky);
  timers.schurBlocksCholesky.stop();
  SchurBlocksCholesky.copyInto(SchurComplementCholesky);

  // TODO: Compute this properly!
  #pragma omp parallel for schedule(dynamic)
  for (int n = 0; n < sdp.polMatrixValues.cols; n++) {
    Real r = sqrt(XInv.diagonalPart[n]*Y.diagonalPart[n]);
    for (int p = 0; p < sdp.polMatrixValues.rows; p++)
      SchurUpdateLowRank.elt(p, n) = r*sdp.polMatrixValues.elt(p, n);
  }

  // cout << "SchurUpdateLowRank[" << i << "] = " << SchurUpdateLowRank << ";\n";

  timers.schurCholeskyUpdate.resume();
  choleskyUpdate(SchurComplementCholesky, SchurUpdateLowRank);
  timers.schurCholeskyUpdate.stop();
}

// v = L^{-1 T}(v - U Q^{-1 T} Q^{-1} U^T v)
//
void partialSchurSolve(BlockDiagonalMatrix &L,
                       const Matrix &U,
                       Matrix &Q,
                       Vector &work,
                       Vector &v) {

  for (int n = 0; n < U.cols; n++) {
    work[n] = 0;
    for (int p = 0; p < U.rows; p++)
      work[n] += U.elt(p, n)*v[p];
  }

  lowerTriangularSolve(Q, work);
  lowerTriangularTransposeSolve(Q, work);

  for (int p = 0; p < U.rows; p++)
    for (int n = 0; n < U.cols; n++)
      v[p] -= U.elt(p, n)*work[n];
      
  blockMatrixLowerTriangularTransposeSolve(L, v);
}

void computeSchurComplementCholesky2(const SDP &sdp,
                                     const BlockDiagonalMatrix &XInv,
                                     const BlockDiagonalMatrix &BilinearPairingsXInv,
                                     const BlockDiagonalMatrix &Y,
                                     const BlockDiagonalMatrix &BilinearPairingsY,
                                     BlockDiagonalMatrix &SchurBlocks,
                                     BlockDiagonalMatrix &SchurBlocksCholesky,
                                     Matrix &SchurUpdateLowRank,
                                     Matrix &P,
                                     Matrix &Q,
                                     Vector &y,
                                     Vector &work) {
  timers.computeSchurBlocks.resume();
  computeSchurBlocks(sdp, BilinearPairingsXInv, BilinearPairingsY, SchurBlocks);
  timers.computeSchurBlocks.stop();

  // temporarily make SchurBlocks full rank
  SchurBlocks.blocks.back().elt(0,0) = 1;

  timers.schurBlocksCholesky.resume();
  choleskyDecomposition(SchurBlocks, SchurBlocksCholesky);
  timers.schurBlocksCholesky.stop();

  // U = SchurBlocksCholesky^{-1} sdp.polMatrixValues
  SchurUpdateLowRank.copyFrom(sdp.polMatrixValues);
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, SchurUpdateLowRank);

  // P = U^T U
  for (int n = 0; n < P.rows; n++) {
    // Parallelize at this loop level because the work is
    // approximately constant size and there are no shared variables;
    #pragma omp parallel for schedule(static)
    for (int m = 0; m <= n; m++) {
      Real tmp = 0;

      for (int p = 0; p < SchurUpdateLowRank.rows; p++)
        tmp += SchurUpdateLowRank.elt(p, n)*SchurUpdateLowRank.elt(p,m);

      P.elt(m, n) = tmp;
      if (m != n)
        P.elt(n, m) = tmp;
    }
  }

  // P = D^{-1} + U^T U
  for (int n = 0; n < P.rows; n++)
    P.elt(n,n) += 1/(XInv.diagonalPart[n]*Y.diagonalPart[n]);

  // P = Q Q^T
  choleskyDecomposition(P, Q);

  // y = u = L^{-1} u = (0, 0, 0, ..., 1)
  fillVector(y, 0);
  y.back() = 1;

  // y = L^{-1 T} (y - U Q^{-1 T} Q^{-1} U^T y)
  partialSchurSolve(SchurBlocksCholesky, SchurUpdateLowRank, Q, work, y);

}

void schurComplementSolve(BlockDiagonalMatrix &SchurBlocksCholesky,
                          const Matrix &SchurUpdateLowRank,
                          Matrix &Q,
                          const Vector &y,
                          Vector &work,
                          Vector &x) {
  // x = L^{-1} x
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, x);

  // x = L^{-1 T} (x - U Q^{-1 T} Q^{-1} U^T x)
  partialSchurSolve(SchurBlocksCholesky, SchurUpdateLowRank, Q, work, x);

  // alpha = (u . x) / (1 - u . y)
  Real alpha = x.back() / (1 - y.back());

  // x = x + alpha y
  for (unsigned int p = 0; p < x.size(); p++)
    x[p] += y[p]*alpha;
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

void SDPSolver::computeSearchDirectionWithRMatrix(const BlockDiagonalMatrix &R,
                                                  const bool isPrimalFeasible) {

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
  blockDiagonalMatrixMultiply(PrimalResidues, Y, Z);
  Z -= R;
  blockMatrixSolveWithInverseCholesky(XInvCholesky, Z);
  Z.symmetrize();

  // dx = schurComplement^-1 r
  computeSchurRHS(sdp, dualResidues, Z, dx);
  vectorSolveWithCholesky(SchurComplementCholesky, dx);
  //schurComplementSolve(SchurBlocksCholesky, SchurUpdateLowRank, schurComplementQ, schurComplementY, schurComplementWork, dx);

  // dX = R_p + sum_p F_p dx_p
  constraintMatrixWeightedSum(sdp, dx, dX);
  dX += PrimalResidues;
  
  // dY = Symmetrize(X^{-1} (R - dX Y))
  blockDiagonalMatrixMultiply(dX, Y, dY);
  dY -= R;
  blockMatrixSolveWithInverseCholesky(XInvCholesky, dY);
  dY.symmetrize();
  dY *= -1;
}

SDPSolverStatus SDPSolver::run(const SDPSolverParameters &parameters,
                               const path outFile,
                               const path checkpointFile) {
  printInfoHeader();
  timers.runSolver.resume();

  for (int iteration = 1; iteration <= parameters.maxIterations; iteration++) {
    inverseCholeskyAndInverse(X, XInvWorkspace, XInvCholesky, XInv);
    inverseCholesky(Y, XInvWorkspace, YInvCholesky);

    timers.bilinearPairings.resume();
    computeBilinearPairings(XInv, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsXInv);
    computeBilinearPairings(Y,    sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsY);
    timers.bilinearPairings.stop();

    // d_k = c_k - Tr(F_k Y)
    computeDualResidues(sdp, Y, BilinearPairingsY, dualResidues);
    // dTilde_N = d_N + E d_B;
    nonBasicShift(E, dualResidues, dualResiduesTilde);

    // y = (x_B - g | x_N)^T
    y = x;
    for (unsigned int i = 0; i < g.size(); i++)
      y[i] -= g[i];
    // PrimalResidues = sum_p F_p x_p - X - F_0
    computePrimalResidues(sdp, y, X, PrimalResidues);

    status.primalError     = PrimalResidues.maxAbsElement();
    status.dualError       = maxAbsVectorElement(dualResidues);
    status.primalObjective = primalObjective(sdp, x);
    status.dualObjective   = dualObjective(sdp, g, dualResidues);
    //status.dualObjective   = dualObjective(sdp, Y);

    const bool isPrimalFeasible = status.isPrimalFeasible(parameters);
    const bool isDualFeasible   = status.isDualFeasible(parameters);
    const bool isOptimal        = status.isOptimal(parameters);
    const bool reductionSwitch  = true;

    if (isPrimalFeasible && isDualFeasible && isOptimal) break;

    timers.schurComplementCholesky.resume();
    computeSchurComplementCholesky(sdp,
                                   XInv, BilinearPairingsXInv,
                                   Y, BilinearPairingsY,
                                   SchurBlocks,
                                   SchurBlocksCholesky,
                                   SchurUpdateLowRank,
                                   SchurComplementCholesky,
                                   iteration);
    // computeSchurComplementCholesky2(sdp,
    //                                 XInv,
    //                                 BilinearPairingsXInv,
    //                                 Y,
    //                                 BilinearPairingsY,
    //                                 SchurBlocks,
    //                                 SchurBlocksCholesky,
    //                                 SchurUpdateLowRank,
    //                                 schurComplementP,
    //                                 schurComplementQ,
    //                                 schurComplementY,
    //                                 schurComplementWork);
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
    Real primalStepLength = stepLength(XInvCholesky, dX, StepMatrixWorkspace,
                                       eigenvaluesWorkspace, QRWorkspace, parameters);
    Real dualStepLength   = stepLength(YInvCholesky, dY, StepMatrixWorkspace,
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
  BlockDiagonalMatrix F(BlockDiagonalMatrix(sdp.objective.size(), sdp.psdMatrixBlockDims()));

  os << "* SDP dense format" << endl;
  os << sdp.affineConstants.size() << " = mDIM" << endl;
  os << F.blocks.size() + 1 << " = nBLOCK" << endl;
  os << "{-" << F.diagonalPart.size();
  for (unsigned int b = 0; b < F.blocks.size(); b++)
    os << ", " << F.blocks[b].rows;
  os << "} = bLOCKsTRUCT" << endl;

  os << sdp.affineConstants << endl;

  F *= 0;
  for (unsigned int n = 0; n < sdp.objective.size(); n++)
    F.diagonalPart[n] = sdp.objective[n];
  os << F << endl;

  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t++) {
      F *= 0;

      for (int n = 0; n < sdp.polMatrixValues.cols; n++)
        F.diagonalPart[n] = sdp.polMatrixValues.elt(t->p, n);

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

  // printSDPBHeader(sdpFile, outFile, checkpointFile, parameters);

  const SDP sdp = readBootstrapSDP(sdpFile);

  // cout << "polMatrixValues = " << sdp.polMatrixValues << ";\n";
  // cout << "bilinearBases = " << sdp.bilinearBases << ";\n";
  
  SDPSolver solver(sdp);
  solver.initialize(parameters);
  // SDPSolverStatus status = solver.run(parameters, outFile, checkpointFile);
  
  // cout << "\nStatus:\n";
  // cout << status << endl;
  // cout << timers << endl;

  // cout << "X = " << solver.X << ";\n";
  // cout << "Y = " << solver.Y << ";\n";
  // cout << "x = " << solver.x << ";\n";
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
  // printSDPDenseFormat(datStream, sdp);
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
  LT = L;
  LT.transposeInplace();

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
  return 0;

  //testBlockCongruence();
  //testBlockDiagonalCholesky();
  //testSDPSolver(argv[1], argv[2]);
  //testCholeskyUpdate();
  //testMinEigenvalue();
}
