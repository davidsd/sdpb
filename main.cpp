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

template <class T>
void printVector(ostream& os, const vector<T> &v) {
  os << "{";
  for (unsigned int i = 0; i < v.size(); i++) {
    os << v[i];
    if (i < v.size() - 1)
      os << ", ";
  }
  os << "}";
}

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
  // A collection of vectors v_k for k=1,...,kmax
  // For example, with vectors of polynomials P_n(x), diagonal Constraints would be
  // {{P_1(x_1),...,P_n(x_1)},...,{P_1(x_kmax),...,P_n(x_kmax)}}
  vector<vector<Real> > diagonalConstraints;
  int row;
  int col;
  vector<int> blocks;

  SDPConstraint(const vector<vector<Real> > &diagonalConstraints,
                const int row,
                const int col,
                const vector<int> &blocks):
    diagonalConstraints(diagonalConstraints),
    row(row),
    col(col),
    blocks(blocks) {}

  static inline SDPConstraint diagonal(const vector<Real> &vec) {
    return SDPConstraint(vector<vector<Real> >(1, vec), 0, 0, vector<int>());
  }

  friend ostream& operator<<(ostream& os, const SDPConstraint& c);

};

ostream& operator<<(ostream& os, const SDPConstraint& c) {
  os << "SDPConstraint(diagonalConstraints = {";
  for (unsigned int i = 0; i < c.diagonalConstraints.size(); i++) {
    printVector(os, c.diagonalConstraints[i]);
    if (i < c.diagonalConstraints.size() - 1)
      os << ", ";
  }
  os << "}, row = " << c.row << ", col = " << c.col << ", blocks = ";
  printVector(os, c.blocks);
  os << ")";

  return os;
}


class SDP {
public:
  // Each algebra basis vector matrix has columns of the form
  // { q_1(x_k), q_2(x_k), ..., q_n(x_k) },
  // with one column for each k = 1,...,kmax.
  //
  // Here, q_i are a bilinear basis for a positive function of some
  // type, in the sense that q_i q_j span the space that the function
  // lives in.  For example, our function might be a polynomial of
  // degree d, in which case we can take the monomial basis q_i = x^i
  // and kmax = d+1 (since we need d+1 values to determine a degree-d
  // polynomial).  The number n of q's needed is d/2+1 if d is even.
  vector<Matrix> algebraBasisVectors;

  // The constraint matrices F_i, in a compressed form.  Each
  // SDPConstraint c represents c.diagonalConstraints.size() different
  // constraints
  vector<SDPConstraint> constraints;

  // The constants c_i on the right-hand side of the affine
  // constraints Tr(F_i Y) = c_i.  We should have
  // affineConstraints.size() = =numConstraints()
  vector<Real> affineConstants;
  vector<Real> objective;

  // Each SDP constraint c represents c.diagonalConstraints.size()
  // different matrices F_i.  This function counts the total number of
  // these matrices.
  int numConstraints() const {
    int result = 0;
    for (unsigned int i = 0; i < constraints.size(); i++) {
      result += constraints[i].diagonalConstraints.size();
    }
    return result;
  }

  // Each algebra basis has n elements, where n is the number of rows
  // in a matrix of algebraBasis values.
  vector<int> blockSizes() const {
    vector<int> sizes;
    for (vector<Matrix>::const_iterator m = algebraBasisVectors.begin();
         m != algebraBasisVectors.end();
         m++)
      sizes.push_back(m->rows);
    return sizes;
  }

  
  // vector<int> () const {
  //   vector<int> sizes;
  //   for (vector<Matrix>::const_iterator m = algebraBasisVectors.begin();
  //        m != algebraBasisVectors.end();
  //        m++)
  //     sizes.push_back(m->rows);
  //   return sizes;
  // }

  friend ostream& operator<<(ostream& os, const SDP& sdp);
};

ostream& operator<<(ostream& os, const SDP& sdp) {
  os << "SDP(algebraBasisVectors = ";
  printVector(os, sdp.algebraBasisVectors);
  os << ", constraints = ";
  printVector(os, sdp.constraints);
  os << ", affineConstants = ";
  printVector(os, sdp.affineConstants);
  os << ", objective = ";
  printVector(os, sdp.objective);
  os << ")";

  return os;
}

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

vector<Real> naturalNumbers(int n) {
  vector<Real> xs(n);
  for (int i = 0; i < n; i++)
    xs[i] = Real(i+1);
  return xs;
}

Matrix monomialAlgebraBasis(int d1, int d, const vector<Real> &xs, bool halfShift) {
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

vector<Real> mapPolynomial(const Polynomial &p, const vector<Real> &xs) {
  vector<Real> ys(xs.size());
  for (unsigned int i = 0; i < xs.size(); i++)
    ys[i] = p(xs[i]);
  return ys;
}

// Given a vector of polynomials {P_1,...,P_n}, and a constant x,
// compute {P_1(x),...,P_n(x)}
//
vector<Real> evalPolVector(const vector<Polynomial> &v, const Real &x) {
  vector<Real> ys(v.size());
  for (unsigned int i = 0; i < v.size(); i++) {
    ys[i] = v[i](x);
  }
  return ys;
}

// Given a vector of polynomials {P_1,...,P_n}, and a vector of
// constants {x_1,...,x_k}, compute the vector of vectors
// {{P_1(x_1),...,P_n(x_1)},...,{P_1(x_k),...,P_n(x_k)}}
//
vector<vector<Real> > evalPolVectors(const vector<Polynomial> &v, const vector<Real> &xs) {
  vector<vector<Real> > ys;
  for (unsigned int k = 0; k < xs.size(); k++) {
    ys.push_back(evalPolVector(v, xs[k]));
  }
  return ys;
}

SDP bootstrapSDP(const vector<Real> &objective,
                 const vector<Real> &normalization,
                 const vector<PolynomialVectorMatrix> &positiveMatrixPols) {
  SDP sdp;
  sdp.objective = objective;
  sdp.constraints.push_back(SDPConstraint::diagonal(normalization));
  sdp.affineConstants.push_back(1);

  for (vector<PolynomialVectorMatrix>::const_iterator m = positiveMatrixPols.begin();
       m != positiveMatrixPols.end();
       m++) {

    int d  = m->degree();
    int d1 = d/2;
    int d2 = (d-1)/2;

    vector<int> blocks;
    vector<Real> xs = naturalNumbers(d+1);

    blocks.push_back(sdp.algebraBasisVectors.size());
    sdp.algebraBasisVectors.push_back(monomialAlgebraBasis(d1, d, xs, false));

    if (d2 >= 0) {
      blocks.push_back(sdp.algebraBasisVectors.size());
      sdp.algebraBasisVectors.push_back(monomialAlgebraBasis(d2, d, xs, true));
    }

    for (int c = 0; c < m->cols; c++) {
      for (int r = 0; r <= c; r++) {
        for (unsigned int k = 0; k < xs.size(); k++)
          sdp.affineConstants.push_back(0);
        sdp.constraints.push_back(SDPConstraint(evalPolVectors(*m->get(r,c), xs), r, c, blocks));
      }
    }

  }

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
  return bootstrapSDP(parseVector(sdpXml->FirstChildElement("objective")->FirstChildElement("vector")),
                      parseVector(sdpXml->FirstChildElement("normalization")->FirstChildElement("vector")),
                      parseMany("polynomialVectorMatrix",
                                parsePolynomialVectorMatrix,
                                sdpXml->FirstChildElement("positiveMatrixPols")));
}

SDP readBootstrapSDP(const char*file) {
  XMLDocument doc;
  doc.LoadFile(file);
  return parseBootstrapSDP(doc.FirstChildElement("sdp"));
}

void testSDP() {
  SDP sdp = readBootstrapSDP("test.sdp");
  cout << sdp << endl;

  const vector<int> blockSizes = sdp.blockSizes();
  BlockDiagonalMatrix y(vector<Real>(sdp.objective.size(), 1), blockSizes);
  BlockDiagonalMatrix x(vector<Real>(sdp.objective.size(), 1), blockSizes);
  y.setIdentity();
  x.setIdentity();
  //  for (int i = 0; i < sdp.algebraBasisVectors)

// void schurComplement(const SDP &sdp,
//                      const BlockDiagonalMatrix &y,
//                      const BlockDiagonalMatrix &xInv,
//                      vector<Matrix> &bilinearPairingsY,
//                      vector<Matrix> &bilinearPairingsXInv,
//                      vector<Matrix> &work,
//                      Matrix &schur) {
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

  //testBlockCongruence();
  //testBlockDiagonalCholesky();
  testSDP();

}
