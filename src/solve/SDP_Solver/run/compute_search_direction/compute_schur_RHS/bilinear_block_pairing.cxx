#include "../../../../SDP_Solver.hxx"

// v^T Z^(blockRow, blockCol) v, where Z^(r,s) is the (r,s)-th dim x
// dim block inside Z.
//
// Input:
// - v        : pointer to the beginning of a vector of length dim
// - dim      : length of the vector v
// - Z        : (k*dim) x (k*dim) matrix, where k > blockRow, blockCol
// - blockRow : integer labeling block of Z
// - blockCol : integer labeling block of Z
// Output: v^T Z^(blockRow, blockCol) v (returned)
//
Real bilinear_block_pairing(const Real *v, const int dim, const Matrix &Z,
                            const int blockRow, const int blockCol)
{
  Real result = 0;

  for(int r = 0; r < dim; r++)
    {
      Real tmp = 0;

      for(int c = 0; c < dim; c++)
        tmp += *(v + c) * Z.elt(blockRow * dim + r, blockCol * dim + c);
      result += *(v + r) * tmp;
    }
  return result;
}
