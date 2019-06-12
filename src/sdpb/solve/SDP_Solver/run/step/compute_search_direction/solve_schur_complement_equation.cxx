#include "lower_triangular_solve.hxx"
#include "../../../../SDP_Solver.hxx"
#include "../../../../lower_triangular_transpose_solve.hxx"

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
void solve_schur_complement_equation(const Block_Diagonal_Matrix &L,
                                     const Block_Matrix &L_inv_B,
                                     const El::DistMatrix<El::BigFloat> &Q,
                                     Block_Vector &dx, Block_Vector &dy)
{
  // Set dx=L^-1 dx
  lower_triangular_solve(L, dx);

  El::DistMatrix<El::BigFloat> dy_dist;
  Zeros(dy_dist, Q.Height(), 1);
  {
    El::Matrix<El::BigFloat> dy_sum;
    Zeros(dy_sum, Q.Height(), 1);

    for(size_t block = 0; block < L_inv_B.blocks.size(); ++block)
      {
        Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
             L_inv_B.blocks[block], dx.blocks[block], El::BigFloat(1),
             dy.blocks[block]);

        // Locally sum contributions to dy
        for(int64_t row = 0; row < dy.blocks[block].LocalHeight(); ++row)
          {
            int64_t global_row(dy.blocks[block].GlobalRow(row));
            for(int64_t column = 0; column < dy.blocks[block].LocalWidth();
                ++column)
              {
                int64_t global_column(dy.blocks[block].GlobalCol(column));
                dy_sum(global_row, global_column)
                  += dy.blocks[block].GetLocal(row, column);
              }
          }
      }

    // Send out updates for dy
    El::BigFloat zero(0);
    for(int64_t row = 0; row < dy_sum.Height(); ++row)
      for(int64_t column = 0; column < dy_sum.Width(); ++column)
        {
          if(dy_sum(row, column) != zero)
            {
              dy_dist.QueueUpdate(row, column, dy_sum(row, column));
            }
        }
  }
  dy_dist.ProcessQueues();

  // Set dy_dist to Q^{-1} dy_dist
  El::cholesky::SolveAfter(El::UpperOrLowerNS::UPPER,
                           El::OrientationNS::NORMAL, Q, dy_dist);
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> dy_local(dy_dist);

  // dx += L_inv_B dy
  for(size_t block = 0; block < L_inv_B.blocks.size(); ++block)
    {
      for(int64_t row = 0; row < dy.blocks[block].LocalHeight(); ++row)
        {
          int64_t global_row(dy.blocks[block].GlobalRow(row));

          for(int64_t column = 0; column < dy.blocks[block].LocalWidth();
              ++column)
            {
              int64_t global_column(dy.blocks[block].GlobalCol(column));
              dy.blocks[block].SetLocal(
                row, column, dy_local.GetLocal(global_row, global_column));
            }
        }
      Gemv(El::OrientationNS::NORMAL, El::BigFloat(1), L_inv_B.blocks[block],
           dy.blocks[block], El::BigFloat(1), dx.blocks[block]);
    }

  // dx = L^{-T} dx
  lower_triangular_transpose_solve(L, dx);
}
