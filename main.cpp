#include <iterator>
#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <assert.h>
#include "types.h"
#include "tinyxml2.h"

using std::vector;
using std::cout;
using std::endl;
using std::ostream;
using std::max;
using std::min;
using std::pair;
using std::ofstream;

using tinyxml2::XMLDocument;
using tinyxml2::XMLElement;

Real Infinity = 1e200;

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

  inline Real get(int r, int c) const {
    return elements[r + c*rows];
  }

  inline void set(int r, int c, const Real &a) {
    elements[r + c*rows] = a;
  }

  inline void addElt(int r, int c, const Real &a) {
    elements[r + c*rows] += a;
  }

  void setZero() {
    fillVector(elements, 0);
  }

  void addDiagonal(const Real &c) {
    assert(rows == cols);

    for (int i = 0; i < rows; i++)
      elements[i * (rows + 1)] += c;
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
        Real tmp = (get(r,c)+get(c,r))/2; 
        set(r, c, tmp);
        set(c, r, tmp);
      }
    }
  }

  void transpose() {
    assert (rows == cols);
    for (int c = 0; c < cols; c++) {
      for (int r = 0; r < c; r++) {
        Real tmp = get(r, c);
        set(r, c, get(c, r));
        set(c, r, tmp);
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
      os << a.get(r,c);
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

// result = A + B
// void matrixAdd(const Matrix &A, const Matrix &B, Matrix &result) {
//   assert(A.cols == B.cols);
//   assert(A.rows == B.rows);
//   assert(A.cols == result.cols);
//   assert(A.rows == result.rows);

//   for (unsigned int i = 0; i < A.elements.size(); i++)
//     result.elements[i] = A.elements[i] + B.elements[i];
// }

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
      result += A.get(r,c)*B.get(r,c);
  result *= 2;

  for (int r = 0; r < A.rows; r++)
    result += A.get(r,r)*B.get(r,r);
  
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
      result += (X.get(r,c) + dX.get(r,c)) * (Y.get(r,c) + dY.get(r,c));
  result *= 2;

  for (int r = 0; r < X.rows; r++)
    result += (X.get(r,r) + dX.get(r,r)) * (Y.get(r,r) + dY.get(r,r));

  return result;
}
 
// Not currently supporting this.  Should probably switch to mpfr...
//
// void matrixMultiplyFirstSym(Matrix &A, Matrix &B, Matrix &result) {
//   assert(A.cols == A.rows);
//   assert(A.cols == B.rows);
//   assert(B.rows == result.rows);
//   assert(B.cols == result.cols);

//   Rsymm("Left", "Upper", B.rows, B.cols, 1,
//         &A.elements[0], A.rows,
//         &B.elements[0], B.rows,
//         0,
//         &result.elements[0], result.rows);
// }

// result = choleskyDecomposition(a) (lower triangular)
// Inputs:
// - a      : dim x dim symmetric matrix
// - result : dim x dim lower-triangular matrix
//
void choleskyDecomposition(Matrix &a, Matrix &result) {
  int dim = a.rows;
  assert(a.cols == dim);
  assert(result.rows == dim);
  assert(result.cols == dim);

  mpackint info;
  Real *resultArray = &result.elements[0];

  Rcopy(dim*dim, &a.elements[0], 1, resultArray, 1);

  // The lower-triangular part of result is now our cholesky matrix
  Rpotrf("Lower", dim, resultArray, dim, &info);
  assert(info == 0);

  // Set the upper-triangular part of the result to zero
  for (int j = 0; j < dim; j++)
    for (int i = 0; i < j; i++)
      result.elements[i + j*dim] = 0;
}

// result = a^-1
// Inputs:
// - a      : dim x dim lower-triangular matrix
// - result : dim x dim lower-triangular matrix
//
void inverseLowerTriangular(Matrix &a, Matrix &result) {
  int dim = a.rows;
  assert(a.cols == dim);
  assert(result.rows == dim);
  assert(result.cols == dim);

  result.setIdentity();
  Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal",
        dim, dim, 1, &a.elements[0], dim, &result.elements[0], dim);
}

// result = choleskyDecomposition(a)^-1
// Inputs:
// - a      : dim x dim symmetric matrix
// - work   : dim x dim matrix
// - result : dim x dim lower-triangular matrix
//
void inverseCholesky(Matrix &a, Matrix &work, Matrix &result) {
  choleskyDecomposition(a, work);
  inverseLowerTriangular(work, result);
}

// b := ACholesky^{-1 T} ACholesky^{-1} b = A^{-1} b
//
// Inputs:
// - ACholesky : dim x dim lower triangular matrix, the Cholesky decomposition of a matrix A
// - b         : vector of length dim (output)
//
void vectorSolveWithCholesky(Matrix &ACholesky, Vector &b) {
  int dim = ACholesky.rows;
  assert(ACholesky.cols == dim);
  assert((int) b.size() == dim);

  Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal",
        dim, 1, 1, &ACholesky.elements[0], dim, &b[0], dim);

  Rtrsm("Left", "Lower", "Transpose", "NonUnitDiagonal",
        dim, 1, 1, &ACholesky.elements[0], dim, &b[0], dim);
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
        tmp += a.get(r, aColOffset + k) * b.get(k, bCol);
      }

      work.set(r, c, tmp);
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
        tmp += b.get(k, bCol) * work.get(workRowOffset + k, c);
      }

      result.set(r, c, tmp);

      // lower triangle is the same as upper triangle
      if (c != r) {
        result.set(c, r, tmp);
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
  b.set(0,0,2);
  b.set(1,0,3);
  b.set(0,1,4);
  b.set(1,1,5);
  b.set(0,2,6);
  b.set(1,2,7);

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

  BlockDiagonalMatrix(int diagonalSize, const vector<int> &blockSizes):
    dim(diagonalSize),
    diagonalPart(Vector(diagonalSize, 0)) {

    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(Matrix(blockSizes[i], blockSizes[i]));
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

    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] += A.blocks[b];
  }

  void operator-=(const BlockDiagonalMatrix &A) {
    addDiagonalPart(A.diagonalPart, -1);

    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] -= A.blocks[b];
  }

  void operator*=(const Real &c) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] *= c;

    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b] *= c;
  }

  void copyFrom(const BlockDiagonalMatrix &A) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] = A.diagonalPart[i];
    
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].copyFrom(A.blocks[b]);
  }

  void symmetrize() {
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

  for (unsigned int b = 0; b < A.blocks.size(); b++)
    lowerTriangularCongruence(A.blocks[b], L.blocks[b]);
}

void inverseCholesky(BlockDiagonalMatrix &A,
                     BlockDiagonalMatrix &work,
                     BlockDiagonalMatrix &AInvCholesky) {
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
    AInvCholesky.diagonalPart[i] = 1/sqrt(A.diagonalPart[i]);

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

  for (unsigned int b = 0; b < A.blocks.size(); b++)
    inverseCholeskyAndInverse(A.blocks[b], work.blocks[b], AInvCholesky.blocks[b], AInv.blocks[b]);
}

void blockMatrixSolveWithInverseCholesky(BlockDiagonalMatrix &AInvCholesky,
                                         BlockDiagonalMatrix &X) {
  for (unsigned int i = 0; i < X.diagonalPart.size(); i++)
    X.diagonalPart[i] *= AInvCholesky.diagonalPart[i] * AInvCholesky.diagonalPart[i];

  for (unsigned int b = 0; b < X.blocks.size(); b++)
    matrixSolveWithInverseCholesky(AInvCholesky.blocks[b], X.blocks[b]);
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
  }

  void addSlack() {
    polMatrixValues.addColumn();
    for (int r = 0; r < polMatrixValues.rows; r++) {
      Real total = 0;
      for (int c = 0; c < polMatrixValues.cols - 1; c++)
        total += polMatrixValues.get(r, c);
      polMatrixValues.set(r, polMatrixValues.cols - 1, -total);
    }
    objective.push_back(-vectorTotal(objective));
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

  const vector<Polynomial> *get(int r, int c) const {
    return &elements[r + c*rows];
  }

  int degree() const {
    int d = 0;
    for (vector<vector<Polynomial> >::const_iterator e = elements.begin(); e != elements.end(); e++)
      for (vector<Polynomial>::const_iterator p = e->begin(); p != e->end(); p++)
        d = max(p->degree(), d);
    return d;
  }

};

Vector naturalNumbers(int n) {
  Vector xs(n);
  for (int i = 0; i < n; i++)
    xs[i] = Real(i+1);
  return xs;
}

Matrix monomialAlgebraBasis(int d1, int d, const Vector &xs, bool halfShift) {
  Matrix basisMatrix(d1+1, d+1);
  for (int k = 0; k < d+1; k++) {
    Real x = xs[k];
    
    Real xToTheN = 1;
    if (halfShift)
      xToTheN = sqrt(x);

    for (int n = 0; n < d1+1; n++) {
      basisMatrix.set(n, k, xToTheN);
      xToTheN *= x;
    }
  }
  return basisMatrix;
}

SDP bootstrapSDP(const Vector &objective,
                 const Vector &normalization,
                 const vector<PolynomialVectorMatrix> &positiveMatrixPols,
                 const Vector &xs) {
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
  sdp.dimensions.push_back(1);
  sdp.degrees.push_back(0);
  numConstraints += 1;

  sdp.polMatrixValues = Matrix(numConstraints, sdp.objective.size());
  sdp.affineConstants = Vector(numConstraints, 0);

  // normalization constraint
  sdp.affineConstants[numConstraints-1] = 1;

  int p = 0;
  for (vector<PolynomialVectorMatrix>::const_iterator m = positiveMatrixPols.begin();
       m != positiveMatrixPols.end();
       m++) {

    int degree = m->degree();
    int delta1 = degree/2;
    int delta2 = (degree-1)/2;

    vector<int> blocks;

    blocks.push_back(sdp.bilinearBases.size());
    sdp.bilinearBases.push_back(monomialAlgebraBasis(delta1, degree, xs, false));

    if (delta2 >= 0) {
      blocks.push_back(sdp.bilinearBases.size());
      sdp.bilinearBases.push_back(monomialAlgebraBasis(delta2, degree, xs, true));
    }

    sdp.blocks.push_back(blocks);

    for (int s = 0; s < m->cols; s++) {
      for (int r = 0; r <= s; r++) {
        for (int k = 0; k <= degree; k++, p++) {
          const Real xk = xs[k];
          for (unsigned int n = 0; n < sdp.objective.size(); n++)
            sdp.polMatrixValues.set(p, n, -(*m->get(r,s))[n](xk));
        }
      }
    }
  }
  assert(p == numConstraints-1);

  // normalization constraint
  for (unsigned int n = 0; n < sdp.objective.size(); n++)
    sdp.polMatrixValues.set(p, n, normalization[n]);
  sdp.blocks.push_back(vector<int>());

  sdp.addSlack();
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
  return parseMany("coord", parseReal, vecXml);
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
                      naturalNumbers(100));
}

SDP readBootstrapSDP(const char*file) {
  XMLDocument doc;
  doc.LoadFile(file);
  return parseBootstrapSDP(doc.FirstChildElement("sdp"));
}

class SolverParameters {
public:
  Real betaStar;
  Real betaBar;
  Real epsilonStar;
  Real epsilonDash;
  Real gammaStar;
  Real alphaMax;
  Real lambdaStar;
  int maxIterations;
  SolverParameters():
    betaStar(0.1),
    betaBar(0.3),
    epsilonStar(1e-30),
    epsilonDash(1e-30),
    gammaStar(0.9),
    alphaMax(100),
    lambdaStar(1e4),
    maxIterations(100) {}
};

class SDPSolver {
public:
  SDP sdp;
  SolverParameters parameters;
  Vector x;
  Vector dx;
  Vector dualResidues;
  Vector XInvYDiag;
  BlockDiagonalMatrix X;
  BlockDiagonalMatrix XInv;
  BlockDiagonalMatrix XInvCholesky;
  BlockDiagonalMatrix Y;
  BlockDiagonalMatrix YInvCholesky;
  BlockDiagonalMatrix Z;
  BlockDiagonalMatrix dX;
  BlockDiagonalMatrix dY;
  BlockDiagonalMatrix Rc;
  BlockDiagonalMatrix PrimalResidues;
  BlockDiagonalMatrix BilinearPairingsXInv;
  BlockDiagonalMatrix BilinearPairingsY;
  Matrix schurComplement;
  Matrix schurComplementCholesky;
  // workspace variables
  BlockDiagonalMatrix XInvWorkspace;
  BlockDiagonalMatrix StepMatrixWorkspace;
  vector<Matrix> bilinearPairingsWorkspace;
  vector<Vector> eigenvaluesWorkspace;
  vector<Vector> QRWorkspace;

  SDPSolver(const SDP &sdp, const SolverParameters &parameters):
    sdp(sdp),
    parameters(parameters),
    x(sdp.numConstraints(), 0),
    dx(x),
    dualResidues(x),
    XInvYDiag(sdp.objective.size(), 0),
    X(sdp.objective.size(), sdp.psdMatrixBlockDims()),
    XInv(X),
    XInvCholesky(X),
    Y(X),
    YInvCholesky(X),
    Z(X),
    dX(X),
    dY(X),
    Rc(X),
    PrimalResidues(X),
    BilinearPairingsXInv(0, sdp.bilinearPairingBlockDims()),
    BilinearPairingsY(BilinearPairingsXInv),
    schurComplement(sdp.numConstraints(), sdp.numConstraints()),
    schurComplementCholesky(schurComplement),
    XInvWorkspace(X),
    StepMatrixWorkspace(X)
  {
    // initialize bilinearPairingsWorkspace, eigenvaluesWorkspace, QRWorkspace 
    for (unsigned int b = 0; b < sdp.bilinearBases.size(); b++) {
      bilinearPairingsWorkspace.push_back(Matrix(X.blocks[b].rows, BilinearPairingsXInv.blocks[b].cols));
      eigenvaluesWorkspace.push_back(Vector(X.blocks[b].rows));
      QRWorkspace.push_back(Vector(3*X.blocks[b].rows - 1));
    }
  }

  void initialize();
  void run();
  void computeSearchDirectionWithRMatrix(const BlockDiagonalMatrix &R, const bool primalFeasible);

};

void computeBilinearPairings(const BlockDiagonalMatrix &A,
                             const vector<Matrix> &bilinearBases,
                             vector<Matrix> &workspace,
                             BlockDiagonalMatrix &result) {
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
        tmp += *(d+n) * V.get(p, n)*V.get(q, n);
      
      result.set(blockRow*V.rows + p, blockCol*V.rows + q, tmp);
      if (p != q)
        result.set(blockRow*V.rows + q, blockCol*V.rows + p, tmp);
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
      tmp += *(v+c) * A.get(blockRow*dim + r, blockCol*dim + c);
    result += *(v+r) * tmp;
  }
  return result;
}

// result = V^T D V, where D=diag(d) is a diagonal matrix
//
// void diagonalCongruence(Real const *d,
//                         const Matrix &V,
//                         const int blockRow,
//                         const int blockCol,
//                         Matrix &result) {
//   int dim = V.rows;
//
//   for (int p = 0; p < V.cols; p++) {
//     for (int q = 0; q <= p; q++) {
//       Real tmp = 0;
//
//       for (int n = 0; n < V.rows; n++)
//         tmp += *(d+n) * V.get(n, p)*V.get(n, q);
//      
//       result.set(blockRow*dim + p, blockCol*dim + q, tmp);
//       if (p != q)
//         result.set(blockRow*dim + q, blockCol*dim + p, tmp);
//     }
//   }
// }

void addSchurBlocks(const SDP &sdp,
                    const BlockDiagonalMatrix &BilinearPairingsXInv,
                    const BlockDiagonalMatrix &BilinearPairingsY,
                    Matrix &schurComplement) {

  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] + 1;

    for (vector<IndexTuple>::const_iterator t1 = sdp.constraintIndices[j].begin();
         t1 != sdp.constraintIndices[j].end();
         t1++) {
      const int p1    = t1->p;
      const int ej_r1 = t1->r * ej;
      const int ej_s1 = t1->s * ej;
      const int k1    = t1->k;

      for (vector<IndexTuple>::const_iterator t2 = sdp.constraintIndices[j].begin();
           t2 != sdp.constraintIndices[j].end() && t2->p <= t1->p;
           t2++) {
        const int p2    = t2->p;
        const int ej_r2 = t2->r * ej;
        const int ej_s2 = t2->s * ej;
        const int k2    = t2->k;

        Real tmp = 0;
        for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
          tmp += (BilinearPairingsXInv.blocks[*b].get(ej_s1 + k1, ej_r2 + k2) * BilinearPairingsY.blocks[*b].get(ej_s2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].get(ej_r1 + k1, ej_r2 + k2) * BilinearPairingsY.blocks[*b].get(ej_s2 + k2, ej_s1 + k1) +
                  BilinearPairingsXInv.blocks[*b].get(ej_s1 + k1, ej_s2 + k2) * BilinearPairingsY.blocks[*b].get(ej_r2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].get(ej_r1 + k1, ej_s2 + k2) * BilinearPairingsY.blocks[*b].get(ej_r2 + k2, ej_s1 + k1))/4;
        }
        schurComplement.addElt(p1, p2, tmp);
        if (p2 != p1)
          schurComplement.addElt(p2, p1, tmp);
      }
    }
  }
}

void computeDualResidues(const SDP &sdp,
                         const BlockDiagonalMatrix &Y,
                         const BlockDiagonalMatrix &BilinearPairingsY,
                         Vector &dualResidues) {
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
        dualResidues[p] -= BilinearPairingsY.blocks[*b].get(ej_r+k, ej_s+k);
        dualResidues[p] -= BilinearPairingsY.blocks[*b].get(ej_s+k, ej_r+k);
      }
      dualResidues[p] /= 2;

      for (int n = 0; n < sdp.polMatrixValues.cols; n++)
        dualResidues[p] -= Y.diagonalPart[n] * sdp.polMatrixValues.get(p, n);

      dualResidues[p] += sdp.affineConstants[p];
    }
  }
}

void constraintMatrixWeightedSum(const SDP &sdp, const Vector x, BlockDiagonalMatrix &result)  {

  for (unsigned int n = 0; n < result.diagonalPart.size(); n++) {
    result.diagonalPart[n] = 0;
    for (unsigned int p = 0; p < x.size(); p++)
      result.diagonalPart[n] += x[p]*sdp.polMatrixValues.get(p, n);
  }

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
    for (unsigned int n = 0; n < Z.diagonalPart.size(); n++)
      r[p] -= sdp.polMatrixValues.get(p, n)*Z.diagonalPart[n];
  }

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

void SDPSolver::initialize() {
  fillVector(x, 0);
  X.setZero();
  X.addDiagonal(parameters.lambdaStar);
  Y.setZero();
  Y.addDiagonal(parameters.lambdaStar);
}

// PrimalResidues = sum_p F_p x_p - X - F_0
//
void computePrimalResidues(const SDP &sdp,
                           const Vector x,
                           const BlockDiagonalMatrix &X,
                           BlockDiagonalMatrix &PrimalResidues) {
  constraintMatrixWeightedSum(sdp, x, PrimalResidues);
  PrimalResidues -= X;
  PrimalResidues.addDiagonalPart(sdp.objective, -1);
}

Real primalObjective(const SDP &sdp, const Vector &x) {
  return dotProduct(sdp.affineConstants, x);
}

Real dualObjective(const SDP &sdp, const BlockDiagonalMatrix &Y) {
  return dotProduct(sdp.objective, Y.diagonalPart);
}

// Implements SDPA's DirectionParameter::MehrotraPredictor
Real predictorCenteringParameter(const SolverParameters &parameters, 
                                 const bool primalFeasible,
                                 const bool dualFeasible,
                                 const bool reductionSwitch) {
  if (primalFeasible && dualFeasible)
    return 0;
  else if (reductionSwitch)
    return parameters.betaBar;
  else
    return 2;
}

// Implements SDPA's DirectionParameter::MehrotraCorrector
Real correctorCenteringParameter(const SolverParameters &parameters,
                                 const BlockDiagonalMatrix &X,
                                 const BlockDiagonalMatrix &dX,
                                 const BlockDiagonalMatrix &Y,
                                 const BlockDiagonalMatrix &dY,
                                 const Real &mu,
                                 const bool primalFeasible,
                                 const bool dualFeasible) {

  Real r = frobeniusProductOfSums(X, dX, Y, dY) / (mu * X.dim);
  Real beta = r < 1 ? r*r : r;

  if (primalFeasible && dualFeasible)
    return min(max(parameters.betaStar, beta), Real(1));
  else
    return max(parameters.betaBar, beta);
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

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
//
Real minEigenvalueViaQR(Matrix &A, Vector &eigenvalues, Vector &workspace) {
  assert(A.rows == A.cols);
  assert((int)eigenvalues.size() == A.rows);
  assert((int)workspace.size() == 3*A.rows - 1);

  mpackint info;
  mpackint workSize = workspace.size();
  Rsyev("NoEigenvectors", "LowerTriangular", A.rows, &A.elements[0], A.rows, &eigenvalues[0], &workspace[0], &workSize, &info);
  assert(info == 0);

  // Eigenvalues are sorted in ascending order
  return eigenvalues[0];
}

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric BlockDiagonalMatrix
// eigenvalues : vector<Vector> of length A.blocks.size()
// workspace   : vector<Vector> of length A.blocks.size()
//
Real minEigenvalueViaQR(BlockDiagonalMatrix &A, vector<Vector> &eigenvalues, vector<Vector> &workspace) {
  assert(A.blocks.size() == eigenvalues.size());
  assert(A.blocks.size() == workspace.size());

  Real lambdaMin = Infinity;
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
    lambdaMin = min(lambdaMin, A.diagonalPart[i]);

  for (unsigned int b = 0; b < A.blocks.size(); b++)
    lambdaMin = min(lambdaMin, minEigenvalueViaQR(A.blocks[b], eigenvalues[b], workspace[b]));

  return lambdaMin;
}

Real stepLength(BlockDiagonalMatrix &XInvCholesky,
                BlockDiagonalMatrix &dX,
                BlockDiagonalMatrix &XInvDX,
                vector<Vector> &eigenvalues,
                vector<Vector> &workspace,
                const SolverParameters &parameters,
                const bool feasible) {

  // XInvDX = L^{-1} dX L^{-1 T}, where X = L L^T
  XInvDX.copyFrom(dX);
  lowerTriangularCongruence(XInvDX, XInvCholesky);
  // cout << "XInvDX = " << XInvDX << endl;

  Real lambda = minEigenvalueViaQR(XInvDX, eigenvalues, workspace);
  cout << "lambda = " << lambda << endl;

  Real alpha = (lambda > -1/parameters.alphaMax) ? parameters.alphaMax : -1/lambda;
  alpha *= parameters.gammaStar;
  if (!feasible)
    alpha = min(alpha, Real(1));
  return alpha;
}

void computeSchurComplementCholesky(const BlockDiagonalMatrix &XInv,
                                    const BlockDiagonalMatrix &BilinearPairingsXInv,
                                    const BlockDiagonalMatrix &Y,
                                    const BlockDiagonalMatrix &BilinearPairingsY,
                                    const SDP &sdp,
                                    Vector &XInvYDiag,
                                    Matrix &schurComplement,
                                    Matrix &schurComplementCholesky) {

  // schurComplement_{pq} = Tr(F_q X^{-1} F_p Y)
  componentProduct(XInv.diagonalPart, Y.diagonalPart, XInvYDiag);

  diagonalCongruence(&XInvYDiag[0], sdp.polMatrixValues, 0, 0, schurComplement);
  addSchurBlocks(sdp, BilinearPairingsXInv, BilinearPairingsY, schurComplement);

  choleskyDecomposition(schurComplement, schurComplementCholesky);

}

void printInfo(Real mu,
               Real primalObj,
               Real dualObj,
               bool primalFeasible,
               bool dualFeasible,
               Real alphaPrimal,
               Real alphaDual) {
  // cout << "---------------------------------------" << endl;
  cout << "mu              = " << mu << endl;
  cout << "primalObjective = " << primalObj << endl;
  cout << "dualObjective   = " << dualObj << endl;
  cout << "primalFeasible  = " << (primalFeasible ? "true" : "false") << endl;
  cout << "dualFeasible    = " << (dualFeasible ? "true" : "false") << endl;
  cout << "alphaPrimal     = " << alphaPrimal << endl;
  cout << "alphaDual       = " << alphaDual << endl;
}

void SDPSolver::computeSearchDirectionWithRMatrix(const BlockDiagonalMatrix &R,
                                                  const bool primalFeasible) {

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R)) if !primalFeasible
  // Z = Symmetrize(X^{-1} (-R))                   if primalFeasible
  if (!primalFeasible)
    blockDiagonalMatrixMultiply(PrimalResidues, Y, Z);
  else
    Z.setZero();
  Z -= R;
  blockMatrixSolveWithInverseCholesky(XInvCholesky, Z);
  Z.symmetrize();

  // dx = schurComplement^-1 r
  computeSchurRHS(sdp, dualResidues, Z, dx);
  vectorSolveWithCholesky(schurComplementCholesky, dx);

  // dX = R_p + sum_p F_p dx_p
  constraintMatrixWeightedSum(sdp, dx, dX);
  // This condition is in the SDPA code but not mentioned in the manual
  if (!primalFeasible)
    dX += PrimalResidues;
  
  // dY = Symmetrize(X^{-1} (R - dX Y))
  blockDiagonalMatrixMultiply(dX, Y, dY);
  dY -= R;
  blockMatrixSolveWithInverseCholesky(XInvCholesky, dY);
  dY.symmetrize();
  dY *= -1;
}

Real dualityGap(Real p, Real d) {
  Real t = (abs(p) + abs(d))/2;
  return abs(p - d) / max(t, Real(1));
}

void SDPSolver::run() {
  for (int iteration = 1; iteration <= parameters.maxIterations; iteration++) {
    cout << "---------------------------------------" << endl;
    cout << "x = " << x << endl;
    cout << "X = " << X << endl;
    cout << "Y = " << Y << endl;

    inverseCholeskyAndInverse(X, XInvWorkspace, XInvCholesky, XInv);
    inverseCholesky(Y, XInvWorkspace, YInvCholesky);

    computeBilinearPairings(XInv, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsXInv);
    computeBilinearPairings(Y,    sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsY);

    // d_k = c_k - Tr(F_k Y)
    computeDualResidues(sdp, Y, BilinearPairingsY, dualResidues);

    // PrimalResidues = sum_p F_p x_p - X - F_0
    computePrimalResidues(sdp, x, X, PrimalResidues);

    bool primalFeasible = PrimalResidues.maxAbsElement()    < parameters.epsilonDash;
    bool dualFeasible   = maxAbsVectorElement(dualResidues) < parameters.epsilonDash;
    Real primalObj      = primalObjective(sdp, x);
    Real dualObj        = dualObjective(sdp, Y);
    bool optimal        = dualityGap(primalObj, dualObj) < parameters.epsilonStar;

    cout << "PrimalResidues.maxAbsElement() = " << PrimalResidues.maxAbsElement() << endl;
    cout << "maxAbsVectorElement(dualResidues) = " << maxAbsVectorElement(dualResidues) << endl;

    bool reductionSwitch = true;

    if (primalFeasible && dualFeasible && optimal)
      return;

    computeSchurComplementCholesky(XInv, BilinearPairingsXInv,
                                   Y,    BilinearPairingsY,
                                   sdp,
                                   XInvYDiag,
                                   schurComplement,
                                   schurComplementCholesky);

    Real mu = frobeniusProductSymmetric(X, Y)/X.dim;

    // Mehrotra predictor solution for (dx, dX, dY)
    Real betaPredictor = predictorCenteringParameter(parameters, primalFeasible, dualFeasible, reductionSwitch);
    cout << "betaPredictor = " << betaPredictor << endl;
    computePredictorRMatrix(betaPredictor, mu, X, Y, Rc);
    computeSearchDirectionWithRMatrix(Rc, primalFeasible);

    cout << "dx_p = " << dx << endl;
    cout << "dX_p = " << dX << endl;
    cout << "dY_p = " << dY << endl;

    // Mehrotra corrector solution for (dx, dX, dY)
    Real betaCorrector = correctorCenteringParameter(parameters, X, dX, Y, dY, mu, primalFeasible, dualFeasible);
    cout << "betaCorrector = " << betaCorrector << endl;
    computeCorrectorRMatrix(betaCorrector, mu, X, dX, Y, dY, Rc);
    computeSearchDirectionWithRMatrix(Rc, primalFeasible);

    cout << "dx_c = " << dx << endl;
    cout << "dX_c = " << dX << endl;
    cout << "dY_c = " << dY << endl;

    // Step length to preserve positive definiteness
    Real alphaPrimal = stepLength(XInvCholesky, dX, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace,
                                  parameters, primalFeasible);
    Real alphaDual   = stepLength(YInvCholesky, dY, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace,
                                  parameters, dualFeasible);

    cout << "alphaPrimal = " << alphaPrimal << endl;
    cout << "alphaDual   = " << alphaDual << endl;

    // Update current point
    vectorScaleMultiplyAdd(alphaPrimal, dx, 1, x);  
    dX *= alphaPrimal;
    X += dX;
    dY *= alphaDual;
    Y += dY;

    if (iteration % 1 == 0)
      printInfo(mu, primalObj, dualObj, primalFeasible, dualFeasible, alphaPrimal, alphaDual);
  }
}

// Compute minimum eigenvalue of L X L^T using the Lanczos method.
// Inputs:
// L : dim x dim Matrix
// X : dim x dim Matrix
// Q : ?   x ?   Matrix
// out     : dim-length Vector
// b       : dim-length Vector
// r       : dim-length Vector
// q       : dim-length Vector
// qold    : dim-length Vector
// w       : dim-length Vector
// tmp     : dim-length Vector
// diagVec : dim-length Vector
// diagVec2: dim-length Vector
// workVec : dim-length Vector
//
Real minEigenvalueViaLanczos(Matrix &L,
                             Matrix &X,
                             Matrix &Q,
                             Vector &out,
                             Vector &b,
                             Vector &r,
                             Vector &q,
                             Vector &qold,
                             Vector &w,
                             Vector &tmp,
                             Vector &diagVec,
                             Vector &diagVec2,
                             Vector &workVec) {
  Real alpha;
  Real value;
  Real min = 1.0e+51;
  Real min_old = 1.0e+52;
  Real min_min= 1.0e+50;
  Real error = 1.0e+10;

  int dim = X.rows;
  int k = 0;
  int kk = 0;
  
  fillVector(diagVec, min_min);
  fillVector(diagVec2, 0);
  fillVector(q, 0);
  fillVector(r, 1);
  
  Real beta = sqrt(Real(dim));  // norm of "r"

  // nakata 2004/12/12
  while (k < dim
         && k < sqrt(Real(dim)) + 10
	 && beta > 1.0e-16
	 && (abs(min-min_old) > (1.0e-5)*abs(min)+(1.0e-8)
	     // && (fabs(min-min_old) > (1.0e-3)*fabs(min)+(1.0e-6)
	     || abs(error*beta) > (1.0e-2)*abs(min)+(1.0e-4) )) {
    cout << "k = " << k << endl;
    cout << "kk = " << kk << endl;

    qold = q;
    value = 1/beta;
    // q = value*r
    vectorScaleMultiplyAdd(value, r, 0, q);

    // w = L X L^T q
    w = q;
    // w = L^T q
    lowerTriangularMatrixTransposeTimesVector(L, w);
    // tmp = X w
    vectorScaleMatrixMultiplyAdd(1, X, w, 0, tmp);
    w = tmp;
    // w = L tmp
    lowerTriangularMatrixTimesVector(L, w);

    alpha = dotProduct(q, w);
    diagVec[k] = alpha;

    // r = w - alpha q - beta qold
    r = w;
    vectorScaleMultiplyAdd(-alpha, q, 1, r);
    vectorScaleMultiplyAdd(-beta, qold, 1, r);

    if ( kk>=sqrt((mpf_class)k) || k==dim-1 || k>sqrt((mpf_class)dim+9) ) {
      kk = 0;
      out = diagVec;
      b = diagVec2;

      out[dim-1] = diagVec[k];
      b[dim-1] = 0;

      mpackint info;
      int kp1 = k+1;
      Rsteqr("I_withEigenvalues", kp1, &out[0], &b[0], &Q.elements[0], Q.rows, &workVec[0], &info);
      
      min_old = min;
      // out have eigen values with ascending order.
      min = out[0];
      error = Q.elements[k];

    }

    value = dotProduct(r,r);
    beta = sqrt(value);
    diagVec2[k] = beta;
    ++k;
    ++kk;

    cout << "beta = " << beta << endl;
    cout << "value = " << value << endl;
  }

  return min - abs(error*beta);
}

void testMinEigenvalue() {
  int dim = 3;

  Matrix L(dim, dim);
  Matrix X(dim, dim);

  L.addDiagonal(1);
  L.set(1,1,2);
  L.set(2,2,3);
  X.addDiagonal(3);
  X.set(1,2,1);
  X.set(2,1,1);
  X.set(0,1,2);
  X.set(1,0,2);

  Matrix Q(dim, dim);
  Vector out(dim);
  Vector b(dim);
  Vector r(dim);
  Vector q(dim);
  Vector qold(dim);
  Vector w(dim);
  Vector tmp(dim);
  Vector diagVec(dim);
  Vector diagVec2(dim);
  Vector workVec(dim);

  Real lambda = minEigenvalueViaLanczos(L, X, Q, out, b, r, q, qold, w, tmp, diagVec, diagVec2, workVec);
  cout << "L = " << L << endl;
  cout << "X = " << X << endl;
  cout << "Q = " << Q << endl;
  cout << "lambda = " << lambda << endl;

  Matrix Y(dim, dim);
  Matrix Work1(dim, dim);
  Matrix Work2(dim, dim);
  Work1 = L;
  Work1.transpose();
  matrixMultiply(X, Work1, Work2);
  matrixMultiply(L, Work2, Y);
  cout << "Y = " << Y << endl;

  Vector Yeigenvalues(dim);
  Vector Yworkspace(3*dim-1);
  cout << "lambdaY = " << minEigenvalueViaQR(Y, Yeigenvalues, Yworkspace) << endl;
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
        F.diagonalPart[n] = sdp.polMatrixValues.get(t->p, n);

      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end();
           b++) {
        const int delta = sdp.bilinearBases[*b].rows;
        const Real *q = &sdp.bilinearBases[*b].elements[(t->k) * delta];

        for (int e = 0; e < delta; e++) {
          for (int f = 0; f < delta; f++) {
            F.blocks[*b].set((t->r)*delta + e, (t->s)*delta + f, (*(q+e)) * (*(q+f)));
          }
        }
        F.blocks[*b].symmetrize();
      }

      os << F << endl;
    }
  }

}

void testSDPSolver(const char *file) {
  const SDP sdp = readBootstrapSDP(file);
  cout << sdp << endl;

  SDPSolver solver(sdp, SolverParameters());
  solver.initialize();
  solver.run();

  cout << "done." << endl;

  // cout << "X = " << solver.X << ";\n";
  // cout << "Y = " << solver.Y << ";\n";
  // cout << "x = " << solver.x << ";\n";
  // cout << "BilinearPairingsXInv = " << solver.BilinearPairingsXInv << endl;
  // cout << "BilinearPairingsY = " << solver.BilinearPairingsY << endl;
  // cout << "schurComplement = " << solver.schurComplement << ";\n";
  // cout << "Rc = " << solver.Rc << ";\n";
  // cout << "dualResidues = " << solver.dualResidues << ";\n";
  // cout << "PrimalResidues = " << solver.PrimalResidues << ";\n";
  // cout << "Z = " << solver.Z << ";\n";
  // cout << "dx = " << solver.dx << ";\n";
  // cout << "dX = " << solver.dX << ";\n";
  // cout << "dY = " << solver.dY << ";\n";

  ofstream datfile;
  datfile.open("sdp.dat");
  printSDPDenseFormat(datfile, sdp);
  datfile.close();
}

int main(int argc, char** argv) {

  mpf_set_default_prec(200);
  cout << "precision = " << mpf_get_default_prec() << endl;
  cout.precision(200);

  //testBlockCongruence();
  //testBlockDiagonalCholesky();
  testSDPSolver(argv[1]);
  //testMinEigenvalue();
  

  return 0;
}
