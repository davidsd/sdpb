#include "../Matrix.hxx"

// y := alpha A x + beta y
void vector_scale_matrix_multiply_add(const Real alpha, const Matrix &A,
                                      const Vector &x, const Real beta,
                                      Vector &y)
{
  assert(A.cols <= x.size());
  assert(A.rows <= y.size());

  for(size_t p = 0; p < A.rows; p++)
    {
      Real tmp = 0;
      for(size_t n = 0; n < A.cols; n++)
        {
          tmp += A.elt(p, n) * x[n];
        }
      y[p] = alpha * tmp + beta * y[p];
    }
}
