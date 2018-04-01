#include "../Matrix.hxx"

// B := L^{-T} B, where L is lower-triangular and B is a matrix
// pointed to by b
void lower_triangular_transpose_solve(Matrix &L, Real *b, int bcols, int ldb)
{
  int dim = L.rows;
  assert(L.cols == dim);

  Rtrsm("Left", "Lower", "Transpose", "NonUnitDiagonal", dim, bcols, 1,
        &L.elements[0], dim, b, ldb);
}

// b := L^{-T} b, where L is lower-triangular
void lower_triangular_transpose_solve(Matrix &L, Vector &b)
{
  assert(static_cast<int>(b.size()) == L.rows);
  lower_triangular_transpose_solve(L, &b[0], 1, b.size());
}
