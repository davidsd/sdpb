#include "../Matrix.hxx"

// Compute the eigenvalues of A, via the QR method
void matrix_eigenvalues(Matrix &A, Vector &workspace, Vector &eigenvalues)
{
  assert(A.rows == A.cols);
  assert(eigenvalues.size() == A.rows);
  assert(workspace.size() == 3 * A.rows - 1);

  Integer info;
  Integer workSize = workspace.size();
  Rsyev("NoEigenvectors", "LowerTriangular", A.rows, A.elements.data(), A.rows,
        eigenvalues.data(), workspace.data(), workSize, &info);
  assert(info == 0);
}
