#include "../Matrix.hxx"

// Set block starting at (bRow, bCol) of B to A^T A
void matrix_square_into_block(Matrix &A, Matrix &B, int bRow, int bCol)
{
  assert(bRow + A.cols <= B.rows);
  assert(bCol + A.cols <= B.cols);

  // This operation is not used within a Block_Diagonal_Matrix, so it is
  // worthwhile to parallelize.  In fact, this function, used in
  // computing TopLeft(Q) = SchurOffDiagonal^T SchurOffDiagonal (see
  // SDPSolver.cpp) is one of the main performance bottlenecks in the
  // solver.

  for(int c = 0; c < A.cols; c++)
    {
      for(int r = 0; r <= c; r++)
        {
          Real tmp = 0;
          for(int p = 0; p < A.rows; p++)
            tmp += A.elt(p, r) * A.elt(p, c);
          B.elt(bRow + r, bCol + c) = tmp;
          if(r != c)
            B.elt(bRow + c, bCol + r) = tmp;
        }
    }
}
