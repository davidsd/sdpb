#include "../Matrix.hxx"

// b := A^{-1} b, where LU and pivots encode the LU decomposition of A
void solve_with_LU_decomposition(Matrix &LU, vector<Integer> &pivots,
                                 Vector &b)
{
  Integer info;
  Rgetrs("NoTranspose", LU.rows, 1, &LU.elements[0], LU.rows, &pivots[0],
         &b[0], b.size(), &info);
  assert(info == 0);
}
