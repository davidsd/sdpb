#include "../Matrix.hxx"

// Compute the eigenvalues of A, via the QR method
void matrix_eigenvalues(Matrix &A, Vector &workspace, Vector &eigenvalues)
{
  assert(A.rows == A.cols);
  assert(static_cast<int>(eigenvalues.size()) == A.rows);
  assert(static_cast<int>(workspace.size()) == 3 * A.rows - 1);

  Integer info;
  Integer workSize = workspace.size();
  Rsyev("NoEigenvectors", "LowerTriangular", A.rows, &A.elements[0], A.rows,
        &eigenvalues[0], &workspace[0], workSize, &info);
  assert(info == 0);
}
