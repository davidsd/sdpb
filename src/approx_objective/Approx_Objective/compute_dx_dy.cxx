#include "sdp_solve/sdp_solve.hxx"

// TODO: Have this be part of sdp_solve.hxx
void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);

void compute_dx_dy(const Block_Info &block_info,
                   const SDP &d_sdp, const Block_Vector &x,
                   const Block_Vector &y,
                   const Block_Diagonal_Matrix &schur_complement_cholesky,
                   const Block_Matrix &schur_off_diagonal,
                   const El::DistMatrix<El::BigFloat> &Q,
                   Block_Vector &dx,
                   Block_Vector &dy)
{
  for(size_t block_index(0); block_index != dx.blocks.size(); ++block_index)
    {
      // dx = -dc + dB.y
      dx.blocks[block_index] = d_sdp.primal_objective_c.blocks[block_index];
      Gemv(El::Orientation::NORMAL, El::BigFloat(1),
           d_sdp.free_var_matrix.blocks[block_index], y.blocks[block_index],
           El::BigFloat(-1), dx.blocks[block_index]);

      // dy = db - x.dB
      El::Zero(dy.blocks[block_index]);
      if(block_info.block_indices[block_index] == 0)
        {
          dy.blocks[block_index] = d_sdp.dual_objective_b;
        }
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               d_sdp.free_var_matrix.blocks[block_index],
               x.blocks[block_index], El::BigFloat(1.0),
               dy.blocks[block_index]);
    }
  // Solve for dx, dy in-place
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal, Q, dx, dy);
}
