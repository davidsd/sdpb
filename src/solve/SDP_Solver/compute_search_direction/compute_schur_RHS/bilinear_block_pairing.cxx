#include "../../../SDP_Solver.hxx"

// v^T A^(blockRow, blockCol) v, where A^(r,s) is the (r,s)-th dim x
// dim block inside A.
//
// Input:
// - v        : pointer to the beginning of a vector of length dim
// - dim      : length of the vector v
// - A        : (k*dim) x (k*dim) matrix, where k > blockRow, blockCol
// - blockRow : integer labeling block of A
// - blockCol : integer labeling block of A
// Output: v^T A^(blockRow, blockCol) v (returned)
//
Real bilinear_block_pairing(const Real *v, const int dim, const Matrix &A,
                            const int blockRow, const int blockCol)
{
  Real result = 0;

  for(int r = 0; r < dim; r++)
    {
      Real tmp = 0;

      for(int c = 0; c < dim; c++)
        tmp += *(v + c) * A.elt(blockRow * dim + r, blockCol * dim + c);
      result += *(v + r) * tmp;
    }
  return result;
}
