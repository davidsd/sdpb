#include "../Matrix.hxx"

// y := alpha A^T x + beta y
void vector_scale_matrix_multiply_transpose_add(Real alpha, Matrix &A,
                                                Vector &x, Real beta,
                                                Vector &y)
{
  assert(A.cols <= static_cast<int>(y.size()));
  assert(A.rows <= static_cast<int>(x.size()));

  for(int n = 0; n < A.cols; n++)
    {
      Real tmp = 0;
      for(int p = 0; p < A.rows; p++)
        tmp += A.elt(p, n) * x[p];
      y[n] = alpha * tmp + beta * y[n];
    }
}
