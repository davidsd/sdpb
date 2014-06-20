#include <iterator>
#include <iostream>
#include <ostream>
#include <vector>
#include <assert.h>
#include "types.h"

using std::vector;
using std::cout;
using std::endl;
using std::ostream;

class Matrix {
 public:
  int rows;
  int cols;
  vector<Real> elements;
  
  Matrix(int rows, int cols):
    rows(rows),
    cols(cols),
    elements(vector<Real>(rows*cols, 0)) {}

  inline Real get(int r, int c) const {
    return elements[r + c*rows];
  }

  inline void set(int r, int c, const Real &a) {
    elements[r + c*rows] = a;
  }

  void setZero() {
    std::fill(elements.begin(), elements.end(), 0);
  }

  void setIdentity() {
    assert(rows == cols);

    setZero();
    for (int i = 0; i < rows; i++)
      elements[i * (rows + 1)] = 1;
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

  // Set the upper-triangular part of the result to zero
  for (int j = 0; j < dim; j++) {
    for (int i = 0; i < j; i++) {
      result.elements[i + j*dim] = 0;
    }
  }
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

// invCholesky = choleskyDecomposition(a)^-1
// inverse = a^-1
// Inputs:
// - a           : dim x dim symmetric matrix
// - work        : dim x dim matrix
// - invCholesky : dim x dim lower-triangular matrix
// - inverse     : dim x dim symmetric matrix
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

// result = b'^T a b', where b' = b \otimes 1
// Inputs:
// - a      : d x d symmetric matrix
// - b      : k x n matrix, where d = numBlocks*k, with numBlocks an integer
// - work   : d x (n*numBlocks) matrix
// - result : (n*numBlocks) x (n*numBlocks) symmetric matrix
//
void blockMatrixCongruence(const Matrix &a, const Matrix &b, Matrix &work, Matrix &result) {
  int numBlocks = a.rows / b.rows;

  assert(result.rows == b.cols * numBlocks);
  assert(result.cols == b.cols * numBlocks);

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
      int bCol       = r % b.cols;
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

class BlockDiagonalMatrix {
public:
  vector<Real> diagonalPart;
  vector<Matrix> blocks;

  BlockDiagonalMatrix(const vector<Real> &diagonalPart, const vector<int> &blockSizes):
    diagonalPart(diagonalPart) {
    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(Matrix(blockSizes[i], blockSizes[i]));
    }
  }

  void setZero() {
    std::fill(diagonalPart.begin(), diagonalPart.end(), 0);

    for (unsigned int b = 0; b < blocks.size(); b++) {
      blocks[b].setZero();
    }
  }

  void setIdentity() {
    std::fill(diagonalPart.begin(), diagonalPart.end(), 1);

    for (unsigned int b = 0; b < blocks.size(); b++) {
      blocks[b].setIdentity();
    }
  }

  friend ostream& operator<<(ostream& os, const BlockDiagonalMatrix& a);

};

ostream& operator<<(ostream& os, const BlockDiagonalMatrix& a) {

  os << "{";
  for (unsigned int i = 0; i < a.diagonalPart.size(); i++) {
    os << a.diagonalPart[i] << ", ";
  }

  for (unsigned int b = 0; b < a.blocks.size(); b++) {
    os << a.blocks[b];
    if (b < a.blocks.size() - 1)
      os << ", ";
  }
  os << "}";
  return os;
}

void inverseCholeskyAndInverse(BlockDiagonalMatrix &a,
                               BlockDiagonalMatrix &work,
                               BlockDiagonalMatrix &invCholesky,
                               BlockDiagonalMatrix &inverse) {

  for (unsigned int i = 0; i < a.diagonalPart.size(); i++) {
    Real d = a.diagonalPart[i];
    invCholesky.diagonalPart[i] = 1/sqrt(d);
    inverse.diagonalPart[i] = 1/d;
  }

  for (unsigned int b = 0; b < a.blocks.size(); b++) {
    inverseCholeskyAndInverse(a.blocks[b],
                              work.blocks[b],
                              invCholesky.blocks[b],
                              inverse.blocks[b]);
  }
}

class SDPConstraint {
public:
  vector<vector<Real> > diagonalConstraints;
  int row;
  int col;
  vector<int> blocks;
};

class SDP {
public:
  vector<Matrix> algebraBasisVectors;
  vector<SDPConstraint> constraints;

  int numConstraints() const {
    int result = 0;
    for (unsigned int i = 0; i < constraints.size(); i++) {
      result += constraints[i].diagonalConstraints.size();
    }
    return result;
  }
};

inline Real quadDotProduct(const vector<Real> &v1,
                           const vector<Real> &v2,
                           const vector<Real> &v3,
                           const vector<Real> &v4) {
  Real result = 0;
  for (unsigned int i = 0; i < v1.size(); i++)
    result += v1[i]*v2[i]*v3[i]*v4[i];
  return result;
}

void schurComplement(const SDP &sdp,
                     const BlockDiagonalMatrix &y,
                     const BlockDiagonalMatrix &xInv,
                     vector<Matrix> &bilinearPairingsY,
                     vector<Matrix> &bilinearPairingsXInv,
                     vector<Matrix> &work,
                     Matrix &schur) {
  assert(schur.rows == sdp.numConstraints());
  assert(schur.rows == schur.cols);

  for (unsigned int b = 0; b < sdp.algebraBasisVectors.size(); b++) {
    blockMatrixCongruence(y.blocks[b],    sdp.algebraBasisVectors[b], work[b], bilinearPairingsY[b]);
    blockMatrixCongruence(xInv.blocks[b], sdp.algebraBasisVectors[b], work[b], bilinearPairingsXInv[b]);
  }

  unsigned int r = 0;
  for (vector<SDPConstraint>::const_iterator f1 = sdp.constraints.begin(); f1 != sdp.constraints.end(); ++f1) {

    int kMax1 = f1->diagonalConstraints.size();
    for (int k1 = 0; k1 < kMax1; k1++, r++) {
      
      int row1 = f1->row;
      int col1 = f1->col;

      unsigned int c = 0;
      for (vector<SDPConstraint>::const_iterator f2 = sdp.constraints.begin(); f2 != sdp.constraints.end(); ++f2) {

        int kMax2 = f2->diagonalConstraints.size();
        for (int k2 = 0; k2 < kMax2; k2++, c++) {

          int row2 = f2->row;
          int col2 = f2->col;

          // TODO: take advantage of the fact that this part of the pairing is symmetric
          Real tmp = quadDotProduct(f1->diagonalConstraints[k1], y.diagonalPart,
                                    f2->diagonalConstraints[k2], xInv.diagonalPart);

          for (unsigned int j1 = 0; j1 < f1->blocks.size(); j1++) {
            int b1 = f1->blocks[j1];

            for (unsigned int j2 = 0; j2 < f2->blocks.size(); j2++) {
              int b2 = f2->blocks[j2];

              if (b1 == b2) {
                int kMax = sdp.algebraBasisVectors[b1].cols;
                assert(kMax == kMax1);
                assert(kMax == kMax2);

                // TODO: make the division by 4 happen only once
                tmp += (bilinearPairingsY[b1].get(row1*kMax+k1,row2*kMax+k2) * bilinearPairingsXInv[b1].get(col2*kMax+k2,col1*kMax+k1) +
                        bilinearPairingsY[b1].get(row1*kMax+k1,col2*kMax+k2) * bilinearPairingsXInv[b1].get(row2*kMax+k2,col1*kMax+k1) +
                        bilinearPairingsY[b1].get(col1*kMax+k1,col2*kMax+k2) * bilinearPairingsXInv[b1].get(row2*kMax+k2,row1*kMax+k1) +
                        bilinearPairingsY[b1].get(col1*kMax+k1,row2*kMax+k2) * bilinearPairingsXInv[b1].get(col2*kMax+k2,row1*kMax+k1))/4;
              }

            }
          }
          
          schur.set(r, c, tmp);

        }
      }

    }

  }

}

void testBlockCongruence() {
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

  blockMatrixCongruence(a, b, work, result);

  cout << a << endl;
  cout << b << endl;
  cout << work << endl;
  cout << result << endl;
  
}

void testBlockDiagonalCholesky() {
  vector<int> sizes;
  sizes.push_back(3);
  sizes.push_back(4);
  vector<Real> diag(2, 0);

  BlockDiagonalMatrix a(diag, sizes);
  a.setIdentity();
  a.diagonalPart[0] = 2;
  a.diagonalPart[1] = 3;
  Real aBlock0[] = {14,3,8,3,10,9,8,9,14};
  a.blocks[0].elements.insert(a.blocks[0].elements.begin(), aBlock0, aBlock0 + 9);

  BlockDiagonalMatrix work(diag, sizes);
  BlockDiagonalMatrix invCholesky(diag, sizes);
  BlockDiagonalMatrix inverse(diag, sizes);

  inverseCholeskyAndInverse(a, work, invCholesky, inverse);

  cout << a << endl;
  cout << invCholesky << endl;
  cout << inverse << endl;
}

int main(int argc, char** argv) {

  testBlockCongruence();
  testBlockDiagonalCholesky();

}
