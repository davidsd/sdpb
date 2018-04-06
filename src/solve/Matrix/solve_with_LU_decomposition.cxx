#include "../Matrix.hxx"

// b := A^{-1} b, where LU and pivots encode the LU decomposition of A
void solve_with_LU_decomposition(Matrix &LU, std::vector<Integer> &pivots,
                                 Vector &b)
{
  Integer info;
  Rgetrs("NoTranspose", LU.rows, 1, LU.elements.data(), LU.rows, pivots.data(),
         b.data(), b.size(), &info);
  assert(info == 0);
}
