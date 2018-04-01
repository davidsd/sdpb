#include "../Matrix.hxx"

// Tr(A B), where A and B are symmetric
Real frobenius_product_symmetric(const Matrix &A, const Matrix &B)
{
  assert(A.rows == B.rows);
  assert(A.cols == B.cols);
  assert(A.rows == A.cols);

  Real result = 0;
  for(int c = 0; c < A.cols; c++)
    for(int r = 0; r < c; r++)
      result += A.elt(r, c) * B.elt(r, c);
  result *= 2;

  for(int r = 0; r < A.rows; r++)
    result += A.elt(r, r) * B.elt(r, r);

  return result;
}
