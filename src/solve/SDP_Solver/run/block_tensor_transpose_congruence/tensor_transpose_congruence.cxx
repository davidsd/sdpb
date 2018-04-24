#include "../../../SDP_Solver.hxx"

// result = base'^T A base', where base' = base \otimes 1, where \otimes
// denotes tensor product and 1 is an mxm identity matrix. Inputs:
// - A      : l*m x l*m symmetric matrix
// - base   : l   x n   matrix
// Workspace:
// - work   : l*m x n*m matrix, intermediate workspace (overwritten)
// Output:
// - result : n*m x n*m symmetric matrix (overwritten)
//
// An explanation of the name: a 'congruence' refers to the action M
// -> A M A^T.  We use 'transpose congruence' to refer to a congruence
// with the transposed matrix M -> A^T M A.  'tensor' here refers to
// the fact that we're performing a congruence with the tensor product
// base \otimes 1.
//
void tensor_transpose_congruence(const Matrix &A, const Matrix &base,
                                 Matrix &work, Matrix &result)
{
  int m = A.rows / base.rows;

  assert(result.rows == base.cols * m);
  assert(result.cols == base.cols * m);

  assert(work.rows == A.rows);
  assert(work.cols == result.cols);

  // work = A base'
  for(size_t c = 0; c < work.cols; c++)
    {
      size_t qCol = c % base.cols;
      size_t aColOffset = (c / base.cols) * base.rows;

      for(size_t r = 0; r < work.rows; r++)
        {
          Real tmp = 0;
          for(size_t k = 0; k < base.rows; k++)
            {
              tmp += A.elt(r, aColOffset + k) * base.elt(k, qCol);
            }

          work.elt(r, c) = tmp;
        }
    }

  // result = base'^T work
  for(size_t c = 0; c < result.cols; c++)
    {
      // since result is symmetric, only compute its upper triangle
      for(size_t r = 0; r <= c; r++)
        {
          size_t qCol = r % base.cols;
          size_t workRowOffset = (r / base.cols) * base.rows;

          Real tmp = 0;
          for(size_t k = 0; k < base.rows; k++)
            {
              tmp += base.elt(k, qCol) * work.elt(workRowOffset + k, c);
            }

          result.elt(r, c) = tmp;

          // lower triangle is the same as upper triangle
          if(c != r)
            {
              result.elt(c, r) = tmp;
            }
        }
    }
}
