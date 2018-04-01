#include "../Matrix.hxx"

// Compute an in-place LU decomposition of A, with pivots, suitable
// for use with 'solveWithLUDecomposition'
void LU_decomposition(Matrix &A, vector<Integer> &pivots)
{
  int dim = A.rows;
  assert(A.cols == dim);

  Integer info;
  Rgetrf(dim, dim, &A.elements[0], dim, &pivots[0], &info);
  assert(info == 0);
}
