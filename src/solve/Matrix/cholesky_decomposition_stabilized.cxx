#include "../Matrix.hxx"

// Compute L (lower triangular) such that A + U U^T = L L^T.
void cholesky_decomposition_stabilized(Matrix &A, Matrix &L,
                                       std::vector<Integer> &stabilize_indices,
                                       std::vector<Real> &stabilize_lambdas,
                                       const double stabilize_threshold)
{
  int dim = A.rows;
  assert(A.cols == dim);
  assert(L.rows == dim);
  assert(L.cols == dim);

  // Set lower-triangular part of L to cholesky decomposition
  L.copy_from(A);
  Integer info;
  RpotrfStabilized("Lower", dim, &L.elements[0], dim, &info, stabilize_indices,
                   stabilize_lambdas, stabilize_threshold);
  assert(info == 0);

  // Set the upper-triangular part of the L to zero
  for(int j = 0; j < dim; j++)
    for(int i = 0; i < j; i++)
      {
        L.elements[i + j * dim] = 0;
      }
}
