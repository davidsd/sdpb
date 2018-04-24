#include "../Matrix.hxx"

// L (lower triangular) such that A = L L^T
void cholesky_decomposition(const Matrix &A, Matrix &L)
{
  size_t dim = A.rows;
  assert(A.cols == dim);
  assert(L.rows == dim);
  assert(L.cols == dim);

  // Set lower-triangular part of L to cholesky decomposition
  L.copy_from(A);
  Integer info;
  Rpotrf("Lower", dim, L.elements.data(), dim, &info);
  assert(info == 0);

  // Set the upper-triangular part of the L to zero
  for(size_t j = 0; j < dim; j++)
    for(size_t i = 0; i < j; i++)
      {
        L.elements[i + j * dim] = 0;
      }
}
