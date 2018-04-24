#include "../Matrix.hxx"

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(Matrix &A, Matrix &L)
{
  size_t dim = A.rows;
  assert(A.cols == dim);
  assert(L.rows == dim);
  assert(L.cols == dim);

  // A := A L^{-T}
  Rtrsm("Right", "Lower", "Transpose", "NonUnitDiagonal", dim, dim, 1,
        L.elements.data(), dim, A.elements.data(), dim);

  // A := L^{-1} A
  Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal", dim, dim, 1,
        L.elements.data(), dim, A.elements.data(), dim);
}
