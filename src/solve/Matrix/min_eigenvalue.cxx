#include "../Matrix.hxx"

// Minimum eigenvalue of A, via the QR method
Real min_eigenvalue(Matrix &A, Vector &workspace, Vector &eigenvalues)
{
  matrix_eigenvalues(A, workspace, eigenvalues);
  return eigenvalues[0];
}
