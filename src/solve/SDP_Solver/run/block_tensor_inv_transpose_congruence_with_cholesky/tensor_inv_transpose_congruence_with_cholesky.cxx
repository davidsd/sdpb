#include "../../../SDP_Solver.hxx"

// Result = Q'^T A^{-1} Q', where Q' = Q \otimes 1, where \otimes
// denotes tensor product and 1 is an mxm idenity matrix.
// Inputs:
// - L      : l*m x l*m cholesky decomposition of A
// - Q      : l   x n   matrix
// Workspace:
// - Work   : l*m x n*m matrix, intermediate workspace (overwritten)
// Output:
// - Result : n*m x n*m symmetric matrix (overwritten)
//
void tensor_inv_transpose_congruence_with_cholesky(const Matrix &L,
                                                   const Matrix &Q,
                                                   Matrix &Work,
                                                   Matrix &Result)
{
  // Work = L^{-1} (Q \otimes 1);
  for(int cw = 0; cw < Work.cols; cw++)
    {
      int mc = cw / Q.cols;

      for(int rw = mc * Q.rows; rw < Work.rows; rw++)
        {
          int mr = rw / Q.rows;

          Real tmp = (mr != mc) ? Real(0) : Q.elt(rw % Q.rows, cw % Q.cols);
          for(int cl = mc * Q.rows; cl < rw; cl++)
            tmp -= L.elt(rw, cl) * Work.elt(cl, cw);

          Work.elt(rw, cw) = tmp / L.elt(rw, rw);
        }
    }

  // Result = Work^T Work
  for(int cr = 0; cr < Result.cols; cr++)
    {
      int mc = cr / Q.cols;

      for(int rr = 0; rr <= cr; rr++)
        {
          int mr = rr / Q.cols;

          Real tmp = 0;
          for(int rw = max(mr, mc) * Q.rows; rw < Work.rows; rw++)
            tmp += Work.elt(rw, cr) * Work.elt(rw, rr);

          Result.elt(rr, cr) = tmp;
          if(rr != cr)
            Result.elt(cr, rr) = tmp;
        }
    }
}
