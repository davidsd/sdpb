#include "../Matrix.hxx"

// X := ACholesky^{-1 T} ACholesky^{-1} X = A^{-1} X
void matrix_solve_with_cholesky(const Matrix &ACholesky, Matrix &X)
{
  int dim = X.rows;
  assert(X.cols == dim);
  assert(ACholesky.rows == dim);
  assert(ACholesky.cols == dim);

  lower_triangular_solve(ACholesky, X.elements.data(), X.cols, X.rows);
  lower_triangular_transpose_solve(ACholesky, X.elements.data(), X.cols,
                                   X.rows);
}
