#include "../../../SDP_Solver.hxx"
#include "../../../lower_triangular_solve.hxx"
#include "../../../lower_triangular_transpose_solve.hxx"

namespace
{
  // X := ACholesky^{-1 T} ACholesky^{-1} X = A^{-1} X
  void vector_solve_with_cholesky(const Matrix &ACholesky, Vector &v)
  {
    assert(ACholesky.rows == v.size());
    assert(ACholesky.cols == v.size());

    lower_triangular_solve(ACholesky, v);
    lower_triangular_transpose_solve(ACholesky, v);
  }
}

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
void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Matrix &schur_off_diagonal, const Matrix &Q, Vector &dx, Vector &dy)
{
  // dx = SchurComplementCholesky^{-1} dx
  lower_triangular_solve(schur_complement_cholesky, dx);

  vector_scale_matrix_multiply_transpose_add(-1, schur_off_diagonal, dx, 1,
                                             dy);
  // dyExtended = Q^{-1} dyExtended
  vector_solve_with_cholesky(Q, dy);

  // dx += SchurOffDiagonal dy
  vector_scale_matrix_multiply_add(1, schur_off_diagonal, dy, 1, dx);

  // dx = SchurComplementCholesky^{-T} dx
  lower_triangular_transpose_solve(schur_complement_cholesky, dx);
}

void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx,
  El::DistMatrix<El::BigFloat> &dy)
{
  // dx = SchurComplementCholesky^{-1} dx
  lower_triangular_solve(schur_complement_cholesky, dx);

  // FIXME: This is terribly serial.
  for(size_t block = 0; block < schur_off_diagonal.blocks.size(); ++block)
    {
      Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
           schur_off_diagonal.blocks[block], dx.blocks[block], El::BigFloat(1),
           dy);
    }

  // dyExtended = Q^{-1} dyExtended
  El::cholesky::SolveAfter(El::UpperOrLowerNS::LOWER,
                           El::OrientationNS::NORMAL, Q, dy);

  // dx += SchurOffDiagonal dy
  // FIXME: This is terribly serial.
  for(size_t block = 0; block < schur_off_diagonal.blocks.size(); ++block)
    {
      Gemv(El::OrientationNS::NORMAL, El::BigFloat(1),
           schur_off_diagonal.blocks[block], dy, El::BigFloat(1),
           dx.blocks[block]);
    }

  // dx = SchurComplementCholesky^{-T} dx
  lower_triangular_transpose_solve(schur_complement_cholesky, dx);
}
