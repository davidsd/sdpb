#ifndef SDP_BOOTSTRAP_DENSEMATRIX_H_
#define SDP_BOOTSTRAP_DENSEMATRIX_H_

#include <vector>
#include "types.h"

using std::vector;

class DenseMatrix {
public:
  int dim;

  // The elements of the matrix, in column-major order, for
  // consistency with lapack's conventions.
  vector<Real> elements;
  
  DenseMatrix(int dim):
    dim(dim),
    elements(vector<Real>(dim*dim, 0)) {}

  inline Real elem(int i, int j) const {
    return elements[i + j*dim];
  }

  // thinking of our matrix M as a block matrix where each block has
  // size v.size() == u.size(), compute the pairing v^T M[ij] u, where
  // M[ij] denotes the i,j-th block of M
  Real blockBilinearForm(int i, int j, const vector<Real> &v, const vector<Real> &u) const {
    int blockDim = v.size();

    Real total = 0;
    for (int k = 0; k < blockDim; k++) {
      Real uk = u[k];
      for (int l = 0; l < blockDim; l++) {
        total += v[l] * uk * elem(i*blockDim+l, j*blockDim + k);
      }
    }

    return total;
  }

  void setToZero() {
    std::fill(elements.begin(), elements.end(), 0);
  }

  void setToIdentity() {
    setToZero();
    for (int i = 0; i < dim; i++)
      elements[i * (dim + 1)] = 1;
  }

  void choleskyDecomposition(DenseMatrix &result) {
    mpackint info;
    Real *resultArray = &result.elements[0];
    Rcopy(dim*dim, &elements[0], 1, resultArray, 1);

    // The lower-triangular part of result is now our cholesky matrix
    Rpotrf("Lower", dim, resultArray, dim, &info);

    // Set the upper-triangular part of the result to zero
    for (int j = 0; j < dim; j++) {
      for (int i = 0; i < j; i++) {
        result.elements[i + j*dim] = 0;
      }
    }
  }

  void inverseLowerTriangular(DenseMatrix &result) {
    result.setToIdentity();
    Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal",
          dim, dim, 1, &elements[0], dim, &result.elements[0], dim);
  }

  void inverseCholesky(DenseMatrix &result, DenseMatrix &workSpace) {
    choleskyDecomposition(workSpace);
    workSpace.inverseLowerTriangular(result);
  }

  void inverseCholeskyAndInverse(DenseMatrix &invCholesky, DenseMatrix &inverse, DenseMatrix &workSpace) {
    inverseCholesky(invCholesky, workSpace);

    inverse.elements = invCholesky.elements;
    Rtrmm("Left", "Lower", "Transpose", "NonUnitDiag", dim, dim, 1,
          &invCholesky.elements[0], dim,
          &inverse.elements[0], dim);
  }

};

#endif  // SDP_BOOTSTRAP_DENSEMATRIX_H_
