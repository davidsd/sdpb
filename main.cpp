#include <iterator>
#include <iostream>
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

using tinyxml2::XMLDocument;
using tinyxml2::XMLElement;

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

Real blockPairing(const SDPConstraint &f1,
                  const SDPConstraint &f2,
                  const unsigned int k1, 
                  const unsigned int k2,
                  const vector<Matrix> &bilinearPairingsY,
                  const vector<Matrix> &bilinearPairingsXInv) {

  int kMax  = f1.diagonalConstraints.size();
  int kMax2 = f2.diagonalConstraints.size();
  assert(kMax == kMax2);

  int i1 = f1.row;
  int j1 = f1.col;
  int i2 = f2.row;
  int j2 = f2.col;

  Real tmp = 0;
  for (vector<int>::const_iterator b1 = f1.blocks.begin(); b1 != f1.blocks.end(); b1++) {
    for (vector<int>::const_iterator b2 = f2.blocks.begin(); b2 != f2.blocks.end(); b2++) {
      if (*b1 == *b2) {
        int b = *b1;
        tmp += (bilinearPairingsY[b].get(i1*kMax+k1,i2*kMax+k2) * bilinearPairingsXInv[b].get(j2*kMax+k2,j1*kMax+k1) +
                bilinearPairingsY[b].get(i1*kMax+k1,j2*kMax+k2) * bilinearPairingsXInv[b].get(i2*kMax+k2,j1*kMax+k1) +
                bilinearPairingsY[b].get(j1*kMax+k1,j2*kMax+k2) * bilinearPairingsXInv[b].get(i2*kMax+k2,i1*kMax+k1) +
                bilinearPairingsY[b].get(j1*kMax+k1,i2*kMax+k2) * bilinearPairingsXInv[b].get(j2*kMax+k2,i1*kMax+k1));
      }
    }
  }
  return tmp/4;
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
  for (vector<SDPConstraint>::const_iterator f1 = sdp.constraints.begin(); f1 != sdp.constraints.end(); f1++) {
    for (unsigned int k1 = 0; k1 < f1->diagonalConstraints.size(); k1++, r++) {

      unsigned int c = 0;
      for (vector<SDPConstraint>::const_iterator f2 = sdp.constraints.begin(); f2 != sdp.constraints.end(); f2++) {
        for (unsigned int k2 = 0; k2 < f2->diagonalConstraints.size(); k2++, c++) {

          schur.set(r, c,
                    blockPairing(*f1, *f2, k1, k2, bilinearPairingsY, bilinearPairingsXInv) +
                    quadDotProduct(f1->diagonalConstraints[k1], y.diagonalPart,
                                   f2->diagonalConstraints[k2], xInv.diagonalPart));

        }
      }
    }
  }
}

class Polynomial {
public:
  vector<Real> coeffs;

  Polynomial(): coeffs(vector<Real>(1, 0)) {}

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

  int degree() {
    int d = 0;
    for (vector<vector<Polynomial> >::const_iterator e = elements.begin(); e != elements.end(); e++)
      for (vector<Polynomial>::const_iterator p = e->begin(); p != e->end(); p++)
        d = max(p->degree(), d);
    return d;
  }

};

SDP bootstrapSDP(const vector<Real> &objective,
                 const vector<Real> &normalization,
                 const vector<PolynomialVectorMatrix> &positiveMatrixPols) {
  return SDP();
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

vector<Real> parseVector(XMLElement *vecXml) {
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
  vector<PolynomialVectorMatrix> positiveMatrixPols =
    parseMany("polynomialVectorMatrix",
              parsePolynomialVectorMatrix,
              sdpXml->FirstChildElement("positiveMatrixPols"));
  vector<Real> objective     = parseVector(sdpXml->FirstChildElement("objective"));
  vector<Real> normalization = parseVector(sdpXml->FirstChildElement("normalization"));
  return bootstrapSDP(objective, normalization, positiveMatrixPols);
}

PolynomialVectorMatrix readPolynomialVectorMatrixFile(const char*file) {
  XMLDocument doc;
  doc.LoadFile(file);
  return parsePolynomialVectorMatrix(doc.FirstChildElement("polynomialVectorMatrix"));
}

void testPolVectorMatrix() {
  PolynomialVectorMatrix m = readPolynomialVectorMatrixFile("test.foo");
  cout << m.elements[0][0] << endl;
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

void testSchurComplement() {
  return;
}

int main(int argc, char** argv) {

  testBlockCongruence();
  testBlockDiagonalCholesky();
  testPolVectorMatrix();

}
