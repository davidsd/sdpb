#include "../Matrix.hxx"

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric Matrices and
// '.' is the Frobenius product.
//
Real frobenius_product_of_sums(const Matrix &X, const Matrix &dX,
                               const Matrix &Y, const Matrix &dY)
{
  Real result = 0;

  for(size_t c = 0; c < X.cols; c++)
    for(size_t r = 0; r < c; r++)
      {
        result += (X.elt(r, c) + dX.elt(r, c)) * (Y.elt(r, c) + dY.elt(r, c));
      }
  result *= 2;

  for(size_t r = 0; r < X.rows; r++)
    {
      result += (X.elt(r, r) + dX.elt(r, r)) * (Y.elt(r, r) + dY.elt(r, r));
    }

  return result;
}
