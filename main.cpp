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

  Matrix(int rows = 0, int cols = 0):
    rows(rows),
    cols(cols),
    elements(vector<Real>(rows*cols, 0)) {}

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
    std::fill(elements.begin(), elements.end(), 0);
  }

  void addIdentity(const Real &c) {
    assert(rows == cols);

    for (int i = 0; i < rows; i++)
      elements[i * (rows + 1)] += c;
  }

  void setIdentity() {
    assert(rows == cols);

    setZero();
    addIdentity(1);
  }

  void scalarMultiply(const Real &c) {
    for (unsigned int i = 0; i < elements.size(); i++)
      elements[i] *= c;
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

void matrixMultiply(Matrix &A, Matrix &B, Matrix &result) {
  assert(A.cols == B.rows);
  assert(A.rows == result.rows);
  assert(B.cols == result.cols);

  Rgemm("N", "N", A.rows, B.cols, A.cols, 1,
        &A.elements[0], A.rows,
        &B.elements[0], B.rows,
        0,
        &result.elements[0], result.rows);
}

Real dotProduct(const vector<Real> &u, const vector<Real> v) {
  Real result = 0;
  for (unsigned int i = 0; i < u.size(); i++)
    result += u[i]*v[i];
  return result;
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
void solveInplaceWithCholesky(Matrix &ACholesky, vector<Real> &b) {
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
void tensorMatrixCongruence(const Matrix &a, const Matrix &b, Matrix &work, Matrix &result) {
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

  tensorMatrixCongruence(a, b, work, result);

  cout << a << endl;
  cout << b << endl;
  cout << work << endl;
  cout << result << endl;
  
}

class BlockDiagonalMatrix {
public:
  int dim;
  vector<Real> diagonalPart;
  vector<Matrix> blocks;

  BlockDiagonalMatrix(int diagonalSize, const vector<int> &blockSizes):
    dim(diagonalSize),
    diagonalPart(vector<Real>(diagonalSize, 0)) {

    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(Matrix(blockSizes[i], blockSizes[i]));
      dim += blockSizes[i];
    }
  }

  void setZero() {
    std::fill(diagonalPart.begin(), diagonalPart.end(), 0);

    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].setZero();
  }

  void addIdentity(const Real &c) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] += c;

    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].addIdentity(c);
  }

  void setIdentity() {
    setZero();
    addIdentity(1);
  }

  void scalarMultiply(const Real &c) {
    for (unsigned int i = 0; i < diagonalPart.size(); i++)
      diagonalPart[i] *= c;

    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].scalarMultiply(c);
  }

  void addDiagonalPart(const vector<Real> &v, const Real &alpha) {
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


  friend ostream& operator<<(ostream& os, const BlockDiagonalMatrix& A);

};

ostream& operator<<(ostream& os, const BlockDiagonalMatrix& A) {

  os << "{";
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++) {
    os << A.diagonalPart[i] << ", ";
  }

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
  Real result = dotProduct(A.diagonalPart, B.diagonalPart);

  for (unsigned int b = 0; b < A.blocks.size(); b++)
    result += frobeniusProductSymmetric(A.blocks[b], B.blocks[b]);

  return result;
}
  

void blockDiagonalMatrixMultiply(BlockDiagonalMatrix &A,
                                 BlockDiagonalMatrix &B,
                                 BlockDiagonalMatrix &result) {
  for (unsigned int i = 0; i < A.diagonalPart.size(); i++)
    result.diagonalPart[i] = A.diagonalPart[i] * B.diagonalPart[i];

  for (unsigned int b = 0; b < A.blocks.size(); b++)
    matrixMultiply(A.blocks[b], B.blocks[b], result.blocks[b]);
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

  for (unsigned int b = 0; b < A.blocks.size(); b++) {
    inverseCholeskyAndInverse(A.blocks[b],
                              work.blocks[b],
                              AInvCholesky.blocks[b],
                              AInv.blocks[b]);
  }
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

class SDP {
public:
  vector<Matrix> bilinearBases;
  int numConstraints;
  int objDimension;
  Matrix polMatrixValues;
  vector<Real> affineConstants;
  vector<Real> objective;
  vector<int> dimensions;
  vector<int> degrees;
  vector<vector<int> > blocks;

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

  friend ostream& operator<<(ostream& os, const SDP& sdp);
};

ostream& operator<<(ostream& os, const SDP& sdp) {
  os << "SDP(bilinearBases = ";
  printVector(os, sdp.bilinearBases);
  os << ", polMatrixValues = " << sdp.polMatrixValues;
  os << ", affineConstants = ";
  printVector(os, sdp.affineConstants);
  os << ", objective = ";
  printVector(os, sdp.objective);
  os << ", dimensions = ";
  printVector(os, sdp.dimensions);
  os << ", degrees = ";
  printVector(os, sdp.degrees);
  os << ", blocks = ";
  os << "{";
  for (vector<vector<int> >::const_iterator b = sdp.blocks.begin();
       b != sdp.blocks.end();
       b++) {
    printVector(os, *b);
    if (b != sdp.blocks.end() - 1)
      os << ", ";
  }
  os << "}";
  os << ")";

  return os;
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

SDP bootstrapSDP(const vector<Real> &objective,
                 const vector<Real> &normalization,
                 const vector<PolynomialVectorMatrix> &positiveMatrixPols,
                 const vector<Real> &xs) {
  SDP sdp;
  sdp.objective = objective;

  sdp.objDimension = objective.size();
  sdp.numConstraints = 0;
  for (vector<PolynomialVectorMatrix>::const_iterator m = positiveMatrixPols.begin();
       m != positiveMatrixPols.end();
       m++) {
    int dimension = m->cols;
    int degree    = m->degree();

    sdp.dimensions.push_back(dimension);
    sdp.degrees.push_back(degree);
    sdp.numConstraints += (degree+1)*dimension*(dimension+1)/2;
  }

  // For the normalization constraint
  sdp.dimensions.push_back(1);
  sdp.degrees.push_back(0);
  sdp.numConstraints += 1;

  sdp.polMatrixValues = Matrix(sdp.numConstraints, sdp.objDimension);
  sdp.affineConstants = vector<Real>(sdp.numConstraints, 0);

  // normalization constraint
  sdp.affineConstants[sdp.numConstraints-1] = 1;

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
          for (int n = 0; n < sdp.objDimension; n++)
            sdp.polMatrixValues.set(p, n, (*m->get(r,s))[n](xk));
        }
      }
    }
  }
  assert(p == sdp.numConstraints-1);

  // normalization constraint
  for (int n = 0; n < sdp.objDimension; n++)
    sdp.polMatrixValues.set(p, n, normalization[n]);
  sdp.blocks.push_back(vector<int>());

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
                                sdpXml->FirstChildElement("positiveMatrixPols")),
                      naturalNumbers(100));
}

SDP readBootstrapSDP(const char*file) {
  XMLDocument doc;
  doc.LoadFile(file);
  return parseBootstrapSDP(doc.FirstChildElement("sdp"));
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

class SolverParameters {
public:
  Real beta;
  SolverParameters(const Real &beta): beta(beta) {}
};

class SDPSolver {
public:
  SDP sdp;
  SolverParameters parameters;
  vector<vector<IndexTuple> > constraintIndexTuples;
  Real mu;
  vector<Real> x;
  vector<Real> dx;
  vector<Real> d;
  vector<Real> XInvYDiag;
  BlockDiagonalMatrix X;
  BlockDiagonalMatrix XInv;
  BlockDiagonalMatrix XInvCholesky;
  BlockDiagonalMatrix Y;
  BlockDiagonalMatrix Z;
  BlockDiagonalMatrix dX;
  BlockDiagonalMatrix dY;
  BlockDiagonalMatrix Rc;
  BlockDiagonalMatrix Rp;
  BlockDiagonalMatrix S;
  BlockDiagonalMatrix T;
  Matrix schurComplement;
  Matrix schurComplementCholesky;
  // workspace variables
  BlockDiagonalMatrix XInvWorkspace;
  vector<Matrix> bilinearPairingsWorkspace;

  SDPSolver(const SDP &sdp, const SolverParameters &parameters):
    sdp(sdp),
    parameters(parameters),
    mu(0),
    x(vector<Real>(sdp.numConstraints, 0)),
    dx(x),
    d(x),
    XInvYDiag(vector<Real>(sdp.objDimension, 0)),
    X(BlockDiagonalMatrix(sdp.objDimension, sdp.psdMatrixBlockDims())),
    XInv(X),
    XInvCholesky(X),
    Y(X),
    Z(X),
    dX(X),
    dY(X),
    Rc(X),
    Rp(X),
    S(BlockDiagonalMatrix(0, sdp.bilinearPairingBlockDims())),
    T(S),
    schurComplement(Matrix(sdp.numConstraints, sdp.numConstraints)),
    schurComplementCholesky(schurComplement),
    XInvWorkspace(X)
  {
    // initialize constraintIndexTuples
    int p = 0;
    for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
      constraintIndexTuples.push_back(vector<IndexTuple>(0));

      for (int s = 0; s < sdp.dimensions[j]; s++) {
        for (int r = 0; r <= s; r++) {
          for (int k = 0; k <= sdp.degrees[j]; k++) {
            constraintIndexTuples[j].push_back(IndexTuple(p, r, s, k));
            p++;
          }
        }
      }
    }

    // initialize bilinearPairingsWorkspace
    for (unsigned int b = 0; b < sdp.bilinearBases.size(); b++)
      bilinearPairingsWorkspace.push_back(Matrix(X.blocks[b].rows, S.blocks[b].cols));
  }

  void computeSearchDirection();

};

void computeBilinearPairings(const BlockDiagonalMatrix &A,
                             const vector<Matrix> &bilinearBases,
                             vector<Matrix> &workspace,
                             BlockDiagonalMatrix &result) {
  for (unsigned int b = 0; b < bilinearBases.size(); b++)
    tensorMatrixCongruence(A.blocks[b], bilinearBases[b], workspace[b], result.blocks[b]);
}
           
// result[i] = u[i] * v[i]
//                
void componentProduct(const vector<Real> &u, const vector<Real> &v, vector<Real> &result) {
  for (unsigned int i = 0; i < u.size(); i++)
    result[i] = u[i] * v[i];
}

// result = V D V^T, where D=diag(d) is a diagonal matrix
//
void diagonalCongruenceTranspose(Real const *d,
                                 const Matrix &V,
                                 const int blockRow,
                                 const int blockCol,
                                 Matrix &result) {
  int dim = V.cols;

  for (int p = 0; p < V.rows; p++) {
    for (int q = 0; q <= p; q++) {
      Real tmp = 0;

      for (int n = 0; n < V.cols; n++)
        tmp += *(d+n) * V.get(p, n)*V.get(q, n);
      
      result.set(blockRow*dim + p, blockCol*dim + q, tmp);
      if (p != q)
        result.set(blockRow*dim + q, blockCol*dim + p, tmp);
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
                    const BlockDiagonalMatrix &S,
                    const BlockDiagonalMatrix &T,
                    const vector<vector<IndexTuple> > &constraintIndexTuples,
                    Matrix &schurComplement) {

  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int mj = sdp.dimensions[j];

    for (vector<IndexTuple>::const_iterator t1 = constraintIndexTuples[j].begin();
         t1 != constraintIndexTuples[j].end();
         t1++) {
      const int p1 = t1->p;
      const int mj_r1 = mj*t1->r;
      const int mj_s1 = mj*t1->s;
      const int k1 = t1->k;

      for (vector<IndexTuple>::const_iterator t2 = constraintIndexTuples[j].begin();
           t2->p <= t1->p && t2 != constraintIndexTuples[j].end();
           t2++) {
        const int p2 = t2->p;
        const int mj_r2 = mj*t2->r;
        const int mj_s2 = mj*t2->s;
        const int k2 = t2->k;

        Real tmp = 0;
        for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
          tmp += (S.blocks[*b].get(mj_s1 + k1, mj_r2 + k2) * T.blocks[*b].get(mj_s2 + k2, mj_r1 + k1) +
                  S.blocks[*b].get(mj_r1 + k1, mj_r2 + k2) * T.blocks[*b].get(mj_s2 + k2, mj_s1 + k1) +
                  S.blocks[*b].get(mj_s1 + k1, mj_s2 + k2) * T.blocks[*b].get(mj_r2 + k2, mj_r1 + k1) +
                  S.blocks[*b].get(mj_r1 + k1, mj_s2 + k2) * T.blocks[*b].get(mj_r2 + k2, mj_s1 + k1))/4;
        }
        schurComplement.addElt(p1, p2, tmp);
        if (p2 != p1)
          schurComplement.addElt(p2, p1, tmp);
      }
    }
  }
}

void computeDVector(const SDP &sdp,
                    const BlockDiagonalMatrix &T,
                    const vector<vector<IndexTuple> > &constraintIndexTuples,
                    vector<Real> &d) {
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int mj = sdp.dimensions[j];

    for (vector<IndexTuple>::const_iterator t = constraintIndexTuples[j].begin();
         t != constraintIndexTuples[j].end();
         t++) {
      const int p = t->p;
      const int mj_r = mj*t->r;
      const int mj_s = mj*t->s;
      const int k = t->k;

      d[p] = 0;
      for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
        d[p] -= T.blocks[*b].get(mj_r+k, mj_s+k);
        d[p] -= T.blocks[*b].get(mj_s+k, mj_r+k);
      }
      d[p] /= 2;
      d[p] += sdp.affineConstants[p];
    }
  }
}

void constraintMatrixWeightedSum(const SDP &sdp, const vector<Real> x, BlockDiagonalMatrix &result)  {

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
          diagonalCongruenceTranspose(&x[p], sdp.bilinearBases[*b], r, s, result.blocks[*b]);
        p += dj + 1;
      }
    }
  }
  assert(p == (int)x.size());

  result.symmetrize();
}

void computeSchurRHS(const SDP &sdp,
                     const vector<vector<IndexTuple> > &constraintIndexTuples,
                     vector<Real> &d,
                     BlockDiagonalMatrix &Z, 
                     vector<Real> &r) {

  for (unsigned int p = 0; p < r.size(); p++) {
    r[p] = -d[p];
    for (unsigned int n = 0; n < Z.diagonalPart.size(); n++)
      r[p] -= sdp.polMatrixValues.get(p, n)*Z.diagonalPart[n];
  }

  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    for (vector<IndexTuple>::const_iterator t = constraintIndexTuples[j].begin();
         t != constraintIndexTuples[j].end();
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

void SDPSolver::computeSearchDirection() {

  X.setIdentity();
  Y.setIdentity();
  mu = frobeniusProductSymmetric(X, Y)/X.dim;

  inverseCholeskyAndInverse(X, XInvWorkspace, XInvCholesky, XInv);

  // schurComplement_{pq} = Tr(F_q X^{-1} F_p Y)
  computeBilinearPairings(XInv, sdp.bilinearBases, bilinearPairingsWorkspace, S);
  computeBilinearPairings(Y,    sdp.bilinearBases, bilinearPairingsWorkspace, T);

  componentProduct(XInv.diagonalPart, Y.diagonalPart, XInvYDiag);

  diagonalCongruenceTranspose(&XInvYDiag[0], sdp.polMatrixValues, 0, 0, schurComplement);
  addSchurBlocks(sdp, S, T, constraintIndexTuples, schurComplement);

  choleskyDecomposition(schurComplement, schurComplementCholesky);

  // Rc = beta mu I - X Y
  blockDiagonalMatrixMultiply(X, Y, Rc);
  Rc.scalarMultiply(-1);
  Rc.addIdentity(parameters.beta*mu);

  // d_k = c_k - Tr(F_k Y)
  computeDVector(sdp, T, constraintIndexTuples, d);

  // Rp = sum_p F_p x_p - X - F_0
  constraintMatrixWeightedSum(sdp, x, Rp);
  Rp -= X;
  Rp.addDiagonalPart(sdp.objective, -1);

  // Z = Symmetrize(X^{-1} (Rp Y - Rc))
  blockDiagonalMatrixMultiply(Rp, Y, Z);
  Z -= Rc;
  blockMatrixSolveWithInverseCholesky(XInvCholesky, Z);
  Z.symmetrize();

  // dx = schurComplement^-1 r
  computeSchurRHS(sdp, constraintIndexTuples, d, Z, dx);
  solveInplaceWithCholesky(schurComplementCholesky, dx);

  // dX = R_p + sum_p F_p dx_p
  constraintMatrixWeightedSum(sdp, dx, dX);
  dX += Rp;
  
  // dY = Symmetrize(X^{-1} (Rc - dX Y))
  blockDiagonalMatrixMultiply(dX, Y, dY);
  dY -= Rc;
  blockMatrixSolveWithInverseCholesky(XInvCholesky, dY);
  dY.symmetrize();
  dY.scalarMultiply(-1);

}

void testSDPSolver() {
  const SDP sdp = readBootstrapSDP("test.sdp");
  cout << sdp << endl;

  SDPSolver solver(sdp, SolverParameters(0.7));
  solver.computeSearchDirection();
  cout << "done." << endl;

  cout << solver.S << endl;
  cout << solver.T << endl;
  cout << "schurComplement: " << solver.schurComplement << endl;
  cout << "Rc = " << solver.Rc << endl;
  cout << "d = ";
  printVector(cout, solver.d);
  cout << endl;
  cout << "Rp = " << solver.Rp << endl;
  cout << "Z = " << solver.Z << endl;
  cout << "dx = ";
  printVector(cout, solver.dx);
  cout << endl;
  cout << "dX = " << solver.dX << endl;
  cout << "dY = " << solver.dY << endl;
}

int main(int argc, char** argv) {

  //testBlockCongruence();
  //testBlockDiagonalCholesky();
  testSDPSolver();

  return 0;
}
