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

  void setZero() {
    std::fill(elements.begin(), elements.end(), 0);
  }

  void addIdentity(const Real &c) {
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
  vector<Real> diagonalPart;
  vector<Matrix> blocks;

  BlockDiagonalMatrix(int diagonalSize, const vector<int> &blockSizes):
    diagonalPart(vector<Real>(diagonalSize, 0)) {
    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(Matrix(blockSizes[i], blockSizes[i]));
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

  void symmetrize() {
    for (unsigned int b = 0; b < blocks.size(); b++)
      blocks[b].symmetrize();
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
  int r;
  int s;
  int k;
  IndexTuple(int r, int s, int k): r(r), s(s), k(k) {}
  IndexTuple() {}
};

class SDPSolver {
public:
  SDP sdp;
  vector<vector<IndexTuple> > constraintIndexTuples;
  vector<Real> r;
  vector<Real> x;
  vector<Real> dx;
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
  // workspace variables
  BlockDiagonalMatrix XInvWorkspace;
  vector<Matrix> bilinearPairingsWorkspace;
  Matrix B;

  SDPSolver(const SDP &sdp):
    sdp(sdp),
    r(vector<Real>(sdp.numConstraints, 0)),
    x(vector<Real>(sdp.numConstraints, 0)),
    dx(vector<Real>(sdp.numConstraints, 0)),
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
    XInvWorkspace(X),
    B(Matrix(sdp.numConstraints, sdp.numConstraints))
  {
    // initialize constraintIndexTuples
    for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
      constraintIndexTuples.push_back(vector<IndexTuple>(0));

      for (int s = 0; s < sdp.dimensions[j]; s++) {
        for (int r = 0; r <= s; r++) {
          for (int k = 0; k <= sdp.degrees[j]; k++) {
            constraintIndexTuples[j].push_back(IndexTuple(r, s, k));
          }
        }
      }
    }

    // initialize bilinearPairingsWorkspace
    for (unsigned int b = 0; b < sdp.bilinearBases.size(); b++)
      bilinearPairingsWorkspace.push_back(Matrix(X.blocks[b].rows, S.blocks[b].cols));
  }

};

void computeBilinearPairings(const BlockDiagonalMatrix &A,
                             const vector<Matrix> &bilinearBases,
                             vector<Matrix> &workspace,
                             BlockDiagonalMatrix &result) {
  for (unsigned int b = 0; b < bilinearBases.size(); b++)
    tensorMatrixCongruence(A.blocks[b], bilinearBases[b], workspace[b], result.blocks[b]);
}
                           
void componentProduct(const vector<Real> &u, const vector<Real> &v, vector<Real> &result) {
  for (unsigned int i = 0; i < u.size(); i++)
    result[i] = u[i] * v[i];
}

void diagonalCongruence(const vector<Real> &d, const Matrix &v, Matrix &result) {
  assert(d.size() == v.cols);
  assert(result.cols == v.rows);
  assert(result.rows == v.rows);

  for (int p = 0; p < v.rows; p++) {
    for (int q = 0; q <= p; q++) {
      Real tmp = 0;

      for (int n = 0; n < v.cols; n++)
        tmp += d[n]*v.get(p, n)*v.get(q, n);
      
      result.set(p, q, tmp);
      if (p != q)
        result.set(q, p, tmp);
    }
  }
}

void testSDPSolver() {
  const SDP sdp = readBootstrapSDP("test.sdp");
  cout << sdp << endl;

  SDPSolver solver(sdp);
  cout << "done." << endl;

  solver.X.setIdentity();
  solver.Y.setIdentity();

  cout << solver.X << endl;
  inverseCholeskyAndInverse(solver.X, solver.XInvWorkspace, solver.XInvCholesky, solver.XInv);

  computeBilinearPairings(solver.XInv, sdp.bilinearBases, solver.bilinearPairingsWorkspace, solver.S);
  computeBilinearPairings(solver.Y,    sdp.bilinearBases, solver.bilinearPairingsWorkspace, solver.T);
  componentProduct(solver.XInv.diagonalPart, solver.Y.diagonalPart, solver.XInvYDiag);

  cout << solver.S << endl;
  cout << solver.T << endl;
}

int main(int argc, char** argv) {

  //testBlockCongruence();
  //testBlockDiagonalCholesky();
  testSDPSolver();

  return 0;
}
