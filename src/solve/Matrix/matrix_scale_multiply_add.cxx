#include "../Matrix.hxx"

// C := alpha*A*B + beta*C
void matrix_scale_multiply_add(Real alpha, Matrix &A, Matrix &B, Real beta,
                               Matrix &C)
{
  assert(A.cols == B.rows);
  assert(A.rows == C.rows);
  assert(B.cols == C.cols);

  Rgemm("N", "N", A.rows, B.cols, A.cols, alpha, &A.elements[0], A.rows,
        &B.elements[0], B.rows, beta, &C.elements[0], C.rows);
}
