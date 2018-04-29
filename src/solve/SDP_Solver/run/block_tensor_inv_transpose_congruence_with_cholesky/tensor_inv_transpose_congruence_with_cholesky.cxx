#include "../../../SDP_Solver.hxx"

// result = bilinear_base^T X^{-1} bilinear_base

void tensor_inv_transpose_congruence_with_cholesky(
  const Matrix &X_cholesky_block, const Matrix &bilinear_base_block,
  Matrix &workspace_block, Matrix &result)
{
  // workspace_block = X_cholesky_block^{-1} (base \otimes 1);
  for(size_t cw = 0; cw < workspace_block.cols; cw++)
    {
      size_t mc = cw / bilinear_base_block.cols;

      for(size_t rw = mc * bilinear_base_block.rows; rw < workspace_block.rows;
          rw++)
        {
          size_t mr = rw / bilinear_base_block.rows;

          Real tmp
            = (mr != mc)
                ? Real(0)
                : bilinear_base_block.elt(rw % bilinear_base_block.rows,
                                          cw % bilinear_base_block.cols);
          for(size_t cl = mc * bilinear_base_block.rows; cl < rw; cl++)
            {
              tmp
                -= X_cholesky_block.elt(rw, cl) * workspace_block.elt(cl, cw);
            }

          workspace_block.elt(rw, cw) = tmp / X_cholesky_block.elt(rw, rw);
        }
    }

  // result = workspace_block^T workspace_block
  for(size_t cr = 0; cr < result.cols; cr++)
    {
      size_t mc = cr / bilinear_base_block.cols;

      for(size_t rr = 0; rr <= cr; rr++)
        {
          size_t mr = rr / bilinear_base_block.cols;

          Real tmp = 0;
          for(size_t rw = max(mr, mc) * bilinear_base_block.rows;
              rw < workspace_block.rows; rw++)
            {
              tmp += workspace_block.elt(rw, cr) * workspace_block.elt(rw, rr);
            }

          result.elt(rr, cr) = tmp;
          if(rr != cr)
            result.elt(cr, rr) = tmp;
        }
    }
}
