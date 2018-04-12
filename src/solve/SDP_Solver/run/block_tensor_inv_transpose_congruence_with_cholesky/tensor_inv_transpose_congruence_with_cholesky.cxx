#include "../../../SDP_Solver.hxx"

// result = (base)'^T A^{-1} (base)', where (base)' = base \otimes 1, where
// \otimes denotes tensor product and 1 is an mxm idenity matrix. Inputs:
// - L      : l*m x l*m cholesky decomposition of A
// - base   : l   x n   matrix
// Workspace:
// - work   : l*m x n*m matrix, intermediate workspace (overwritten)
// Output:
// - result : n*m x n*m symmetric matrix (overwritten)
//
void tensor_inv_transpose_congruence_with_cholesky(const Matrix &L,
                                                   const Matrix &base,
                                                   Matrix &work,
                                                   Matrix &result)
{
  // work = L^{-1} (base \otimes 1);
  for(int cw = 0; cw < work.cols; cw++)
    {
      int mc = cw / base.cols;

      for(int rw = mc * base.rows; rw < work.rows; rw++)
        {
          int mr = rw / base.rows;

          Real tmp
            = (mr != mc) ? Real(0) : base.elt(rw % base.rows, cw % base.cols);
          for(int cl = mc * base.rows; cl < rw; cl++)
            tmp -= L.elt(rw, cl) * work.elt(cl, cw);

          work.elt(rw, cw) = tmp / L.elt(rw, rw);
        }
    }

  // result = work^T work
  for(int cr = 0; cr < result.cols; cr++)
    {
      int mc = cr / base.cols;

      for(int rr = 0; rr <= cr; rr++)
        {
          int mr = rr / base.cols;

          Real tmp = 0;
          for(int rw = max(mr, mc) * base.rows; rw < work.rows; rw++)
            tmp += work.elt(rw, cr) * work.elt(rw, rr);

          result.elt(rr, cr) = tmp;
          if(rr != cr)
            result.elt(cr, rr) = tmp;
        }
    }
}
