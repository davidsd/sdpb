#include "../Matrix.hxx"

// y := alpha A x + beta y
void vector_scale_matrix_multiply_add(Real alpha, Matrix &A, Vector &x,
                                      Real beta, Vector &y)
{
  assert(A.cols <= static_cast<int>(x.size()));
  assert(A.rows <= static_cast<int>(y.size()));

  for(int p = 0; p < A.rows; p++)
    {
      Real tmp = 0;
      for(int n = 0; n < A.cols; n++)
        {
          tmp += A.elt(p, n) * x[n];
        }
      y[p] = alpha * tmp + beta * y[p];
    }
}
