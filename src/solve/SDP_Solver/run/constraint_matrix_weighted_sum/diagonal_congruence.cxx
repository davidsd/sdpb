#include "../../../SDP_Solver.hxx"

// Result^(blockRow,blockCol) = V D V^T, where D=diag(d) is a diagonal
// matrix.
//
// Here, we view Result as a matrix on the tensor product
// R^V.rows \otimes R^k.  Result^(blockRow,blockCol) refers to the
// block submatrix of size V.rows x V.rows at position (blockRow,
// blockCol) in the second tensor factor.
//
// Inputs:
// - d        : pointer to beginning of a length-V.cols vector
// - V        : V.rows x V.cols Matrix
// - blockRow : integer < k
// - blockCol : integer < k
// Output:
// - Result   : (k*V.rows) x (k*V.rows) square Matrix (overwritten)
//
void diagonal_congruence(const Real *d, const Matrix &V, const int blockRow,
                         const int blockCol, Matrix &Result)
{
  for(int p = 0; p < V.rows; p++)
    {
      for(int q = 0; q <= p; q++)
        {
          Real tmp = 0;

          for(int n = 0; n < V.cols; n++)
            {
              tmp += *(d + n) * V.elt(p, n) * V.elt(q, n);
            }

          Result.elt(blockRow * V.rows + p, blockCol * V.rows + q) = tmp;
          if(p != q)
            {
              Result.elt(blockRow * V.rows + q, blockCol * V.rows + p) = tmp;
            }
        }
    }
}
