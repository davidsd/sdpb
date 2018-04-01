#include "../SDP_Solver.hxx"

// Solve the Schur complement equation for dx, dy.
//
// - As inputs, dx and dy are the residues r_x and r_y on the
//   right-hand side of the Schur complement equation.
// - As outputs, dx and dy are overwritten with the solutions of the
//   Schur complement equation.
//
// The equation is solved using the block-decomposition described in
// the manual.
//
void SDP_Solver::solve_schur_complement_equation(Vector &dx, Vector &dy)
{
  // dx = SchurComplementCholesky^{-1} dx
  block_matrix_lower_triangular_solve(schur_complement_cholesky, dx);

  // extend dy by additional coordinates z needed for stabilizing
  // dyExtended = (dy, z)
  dy_extended.resize(Q.rows);

  // dy = r_y - SchurOffDiagonal^T dx
  for(size_t n = 0; n < dy.size(); n++)
    {
      dy_extended[n] = dy[n];
    }
  vector_scale_matrix_multiply_transpose_add(-1, schur_off_diagonal, dx, 1,
                                             dy_extended);

  // dyExtended = Q^{-1} dyExtended
  solve_with_LU_decomposition(Q, Q_pivots, dy_extended);

  // dx += SchurOffDiagonal dy
  vector_scale_matrix_multiply_add(1, schur_off_diagonal, dy_extended, 1, dx);

  // dx = SchurComplementCholesky^{-T} dx
  block_matrix_lower_triangular_transpose_solve(schur_complement_cholesky, dx);

  // dy = first few entries of dyExtended
  for(unsigned int n = 0; n < dy.size(); n++)
    dy[n] = dy_extended[n];
}
