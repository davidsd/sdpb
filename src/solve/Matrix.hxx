//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Vector.hxx"

// A matrix M with Real entries
class Matrix
{
public:
  int rows;
  int cols;
  // Elements of M in row-major order:
  //
  //   elements = { M_{0,0}, ..., M_{cols-1,0},
  //                M_{0,1}, ..., M_{cols-1,1},
  //                ...
  //                M_{0,rows-1}, ..., M_{cols-1,rows-1} }
  //
  Vector elements;

  Matrix(int rows = 0, int cols = 0)
      : rows(rows), cols(cols), elements(Vector(rows * cols, 0))
  {}

  Matrix(int rows, int cols, Real *array)
      : rows(rows), cols(cols), elements(Vector(array, array + rows * cols))
  {}

  inline const Real &elt(const int r, const int c) const
  {
    return elements[r + c * rows];
  }

  inline Real &elt(const int r, const int c) { return elements[r + c * rows]; }

  // M := 0
  void set_zero() { std::fill(elements.begin(), elements.end(), 0); }

  // M += c*I, where I is the identity and c is a constant
  void add_diagonal(const Real &c)
  {
    // ensure M is square
    assert(rows == cols);

    for(int i = 0; i < rows; i++)
      {
        elt(i, i) += c;
      }
  }

  void resize(int r, int c)
  {
    elements.resize(r * c);
    rows = r;
    cols = c;
  }

  void symmetrize()
  {
    assert(rows == cols);

    for(int r = 0; r < rows; r++)
      {
        for(int c = 0; c < r; c++)
          {
            Real tmp = (elt(r, c) + elt(c, r)) / 2;
            elt(r, c) = tmp;
            elt(c, r) = tmp;
          }
      }
  }

  // M := A

  // The difference from the default assignment operator is that this
  // never reallocates.  I am not sure it is really needed.
  void copy_from(const Matrix &A)
  {
    assert(rows == A.rows);
    assert(cols == A.cols);

    std::copy(A.elements.begin(), A.elements.end(), elements.begin());
  }

  // M := M + A
  void operator+=(const Matrix &A)
  {
    for(size_t i = 0; i < elements.size(); i++)
      {
        elements[i] += A.elements[i];
      }
  }

  // M := M - A
  void operator-=(const Matrix &A)
  {
    for(size_t i = 0; i < elements.size(); i++)
      {
        elements[i] -= A.elements[i];
      }
  }

  // M := c*M, where c is a constant
  void operator*=(const Real &c)
  {
    for(auto &element : elements)
      {
        element *= c;
      }
  }

  Matrix operator-(const Matrix &other) const
  {
    Matrix out(rows, cols);
    std::transform(elements.begin(), elements.end(), other.elements.begin(),
                   out.elements.begin(), std::minus<Real>());
    return out;
  }

  // The maximum absolute value of the elemnts of M
  Real max_abs() const { return max_abs_vector(elements); }

  friend std::ostream &operator<<(std::ostream &os, const Matrix &a);
};

// Most of the routines below are helpfully-named wrappers for
// functions in MBLAS or MLAPACK.  See Matrix.h for a more detailed
// description of the various input/output parameters.
//
// For a list of MBLAS routines with documentation, see
// http://mplapack.sourceforge.net/mblas_routines.html
//
// For a list of MLAPACK routines with documentation, see
// http://mplapack.sourceforge.net/mlapack_routines.html
//
// We have chosen not to parallelize operations that are used in
// Block_Diagonal_Matrix, since there parallelism can be achieved by
// parallelizing loops over blocks.

// C := alpha*A*B + beta*C
void matrix_scale_multiply_add(Real alpha, Matrix &A, Matrix &B, Real beta,
                               Matrix &C);

// Set block starting at (bRow, bCol) of B to A^T A
void matrix_square_into_block(Matrix &A, Matrix &B, int bRow, int bCol);

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(Matrix &A, Matrix &L);

// y := alpha A x + beta y
void vector_scale_matrix_multiply_add(const Real alpha, const Matrix &A,
                                      const Vector &x, const Real beta,
                                      Vector &y);

// y := alpha A^T x + beta y
void vector_scale_matrix_multiply_transpose_add(const Real alpha,
                                                const Matrix &A,
                                                const Vector &x,
                                                const Real beta, Vector &y);

// Frobenius product Tr(A^T B) where A and B are symmetric matrices
Real frobenius_product_symmetric(const Matrix &A, const Matrix &B);

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric Matrices and
// '.' is the Frobenius product.
Real frobenius_product_of_sums(const Matrix &X, const Matrix &dX,
                               const Matrix &Y, const Matrix &dY);

// Eigenvalues of A, via the QR method
// Inputs:
// A           : n x n Matrix (will be overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
void matrix_eigenvalues(Matrix &A, Vector &workspace, Vector &eigenvalues);

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : n x n Matrix (overwritten)
// eigenvalues : Vector of length n
// workspace   : Vector of lenfth 3*n-1 (temporary workspace)
Real min_eigenvalue(Matrix &A, Vector &workspace, Vector &eigenvalues);

// Compute an in-place LU decomposition of A, with pivots, suitable
// for use with 'solveWithLUDecomposition'
void LU_decomposition(Matrix &A, std::vector<Integer> &pivots);

// b := A^{-1} b, where LU and pivots encode the LU decomposition of A
void solve_with_LU_decomposition(const Matrix &LU,
                                 const std::vector<Integer> &pivots,
                                 Vector &b);

// L (lower triangular) such that A = L L^T
// Input:
// - A : dim x dim symmetric matrix
// Output:
// - L : dim x dim lower-triangular matrix (overwritten)
void cholesky_decomposition(Matrix &A, Matrix &L);

// B := L^{-1} B, where L is lower-triangular and B is a matrix
// pointed to by b
//
// Input:
// - L
// - b, a pointer to the top-left element of B
// - bcols, the number of columns of B
// - ldb, the distance in pointer-space between the first elements of
//   each row of B
// Output:
// - elements of B, which are values pointed to by b are overwritten
//
// (The use of pointers and ldb is to allow using this function for
// submatrices of a larger matrix.)
void lower_triangular_solve(const Matrix &L, Real *b, int bcols, int ldb);

// b := L^{-1} b, where L is lower-triangular
void lower_triangular_solve(const Matrix &L, Vector &b);

// B := L^{-T} B, where L is lower-triangular and B is a matrix
// pointed to by b
//
// Input:
// - L
// - b, a pointer to the top-left element of B
// - bcols, the number of columns of B
// - ldb, the distance in pointer-space between the first elements of
//   each row of B
// Output:
// - elements of B, which are values pointed to by b are overwritten
//
void lower_triangular_transpose_solve(const Matrix &L, Real *b, int bcols,
                                      int ldb);

// b := L^{-T} b, where L is lower-triangular
void lower_triangular_transpose_solve(const Matrix &L, Vector &b);

// X := ACholesky^{-1 T} ACholesky^{-1} X = A^{-1} X
// Inputs:
// - ACholesky : dim x dim lower triangular matrix (constant)
// Output:
// - X         : dim x dim matrix (overwritten)
//
void matrix_solve_with_cholesky(const Matrix &ACholesky, Matrix &X);
