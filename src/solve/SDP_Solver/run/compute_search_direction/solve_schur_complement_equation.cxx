#include "../../../SDP_Solver.hxx"
#include "../../../lower_triangular_solve.hxx"
#include "../../../lower_triangular_transpose_solve.hxx"

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
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy)
{
  // dx = SchurComplementCholesky^{-1} dx
  lower_triangular_solve(schur_complement_cholesky, dx);

  El::DistMatrix<El::BigFloat> dy_dist;
  Zeros(dy_dist, Q.Height(), 1);

  dy_dist.Reserve(schur_off_diagonal.blocks.size() * dy_dist.Height());
  for(size_t block = 0; block < schur_off_diagonal.blocks.size(); ++block)
    {
      Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
           schur_off_diagonal.blocks[block], dx.blocks[block], El::BigFloat(1),
           dy.blocks[block]);

      for(size_t row = 0; row < dy.blocks[block].Height(); ++row)
        for(size_t column = 0; column < dy.blocks[block].Width(); ++column)
          {
            // FIXME: This assumes that there is only one processor
            // per grid.  When fixing this, remember to also fix the
            // Reserve() call.
            dy_dist.QueueUpdate(row, column,
                                dy.blocks[block].GetLocal(row, column));
          }
    }
  dy_dist.ProcessQueues();

  // dyExtended = Q^{-1} dyExtended
  El::cholesky::SolveAfter(El::UpperOrLowerNS::LOWER,
                           El::OrientationNS::NORMAL, Q, dy_dist);
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> dy_local(dy_dist);

  // dx += SchurOffDiagonal dy
  for(size_t block = 0; block < schur_off_diagonal.blocks.size(); ++block)
    {
      for(size_t row = 0; row < dy.blocks[block].Height(); ++row)
        for(size_t column = 0; column < dy.blocks[block].Width(); ++column)
          {
            dy.blocks[block].SetLocal(row, column,
                                      dy_local.GetLocal(row, column));
          }

      Gemv(El::OrientationNS::NORMAL, El::BigFloat(1),
           schur_off_diagonal.blocks[block], dy.blocks[block], El::BigFloat(1),
           dx.blocks[block]);
    }

  // dx = SchurComplementCholesky^{-T} dx
  lower_triangular_transpose_solve(schur_complement_cholesky, dx);
}
