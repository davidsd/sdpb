#include "../Matrix.hxx"

// C := alpha*A*B + beta*C
void matrix_scale_multiply_add(const Real &alpha, Matrix &A,
                               Matrix &B, const Real &beta, Matrix &C)
{
  assert(A.cols == B.rows);
  assert(A.rows == C.rows);
  assert(B.cols == C.cols);

  Rgemm("N", "N", A.rows, B.cols, A.cols, alpha, A.elements.data(), A.rows,
        B.elements.data(), B.rows, beta, C.elements.data(), C.rows);
}
