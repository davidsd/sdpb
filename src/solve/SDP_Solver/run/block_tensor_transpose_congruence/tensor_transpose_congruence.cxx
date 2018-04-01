#include "../../../SDP_Solver.hxx"

// Result = Q'^T A Q', where Q' = Q \otimes 1, where \otimes denotes
// tensor product and 1 is an mxm identity matrix.
// Inputs:
// - A      : l*m x l*m symmetric matrix
// - Q      : l   x n   matrix
// Workspace:
// - Work   : l*m x n*m matrix, intermediate workspace (overwritten)
// Output:
// - Result : n*m x n*m symmetric matrix (overwritten)
//
// An explanation of the name: a 'congruence' refers to the action M
// -> A M A^T.  We use 'transpose congruence' to refer to a congruence
// with the transposed matrix M -> A^T M A.  'tensor' here refers to
// the fact that we're performing a congruence with the tensor product
// Q \otimes 1.
//
void tensor_transpose_congruence(const Matrix &A, const Matrix &Q,
                                 Matrix &Work, Matrix &Result)
{
  int m = A.rows / Q.rows;

  assert(Result.rows == Q.cols * m);
  assert(Result.cols == Q.cols * m);

  assert(Work.rows == A.rows);
  assert(Work.cols == Result.cols);

  // Work = A Q'
  for(int c = 0; c < Work.cols; c++)
    {
      int qCol = c % Q.cols;
      int aColOffset = (c / Q.cols) * Q.rows;

      for(int r = 0; r < Work.rows; r++)
        {
          Real tmp = 0;
          for(int k = 0; k < Q.rows; k++)
            {
              tmp += A.elt(r, aColOffset + k) * Q.elt(k, qCol);
            }

          Work.elt(r, c) = tmp;
        }
    }

  // Result = Q'^T Work
  for(int c = 0; c < Result.cols; c++)
    {
      // since Result is symmetric, only compute its upper triangle
      for(int r = 0; r <= c; r++)
        {
          int qCol = r % Q.cols;
          int workRowOffset = (r / Q.cols) * Q.rows;

          Real tmp = 0;
          for(int k = 0; k < Q.rows; k++)
            {
              tmp += Q.elt(k, qCol) * Work.elt(workRowOffset + k, c);
            }

          Result.elt(r, c) = tmp;

          // lower triangle is the same as upper triangle
          if(c != r)
            {
              Result.elt(c, r) = tmp;
            }
        }
    }
}
