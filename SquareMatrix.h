#ifndef SDP_BOOTSTRAP_SQUAREMATRIX_H_
#define SDP_BOOTSTRAP_SQUAREMATRIX_H_

#include <vector>
#include <assert.h>
#include <ostream>
#include "types.h"

using std::vector;
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

void blockMatrixCongruence(const Matrix &a, const Matrix &b, Matrix &work, Matrix &result);

void choleskyDecomposition(Matrix &a, Matrix &result);

void inverseLowerTriangular(Matrix &a, Matrix &result);

void inverseCholesky(Matrix &a, Matrix &work, Matrix &result);

void inverseCholeskyAndInverse(Matrix &a, Matrix &work, Matrix &invCholesky, Matrix &inverse);

#endif  // SDP_BOOTSTRAP_SQUAREMATRIX_H_
