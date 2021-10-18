#include "sdp_solve/sdp_solve.hxx"

void scale_block_vector(Block_Vector &A, const El::BigFloat &alpha)
{
  for(size_t block = 0; block != A.blocks.size(); ++block)
    {
      A.blocks[block] *= alpha;
    }
}

void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);

//Given delta_p(sdp) , compute the (delta_x, delta_y) = H_xx^-1 H_xp delta_p
//as shown in Eq(15).
//Return: void.
//Update: hess_xp = H_xp = (RHS(p1)/alpha, RHS(p2)/alpha, ... ), stored to compute the second term on the LHS of Eq(13)
//        delta_x_y = - H^-1_xx H_xp = - H^-1_xx hess_xp
void mixed_hess(const Block_Info &block_info, const SDP &d_sdp,
                const Block_Vector &x, const Block_Vector &y,
                const Block_Diagonal_Matrix &schur_complement_cholesky,
                const Block_Matrix &schur_off_diagonal,
                const El::DistMatrix<El::BigFloat> &Q,
                const El::BigFloat &alpha,
                std::vector<std::pair<Block_Vector, Block_Vector>> &hess_xp,
                std::vector<std::pair<Block_Vector, Block_Vector>> &delta_x_y)
{
  Block_Vector dx(x), dy(y);

  for(size_t block_index(0); block_index != dx.blocks.size(); ++block_index)
    {
      // dx = -dc + dB.y
      dx.blocks[block_index] = d_sdp.primal_objective_c.blocks[block_index];
      Gemv(El::Orientation::NORMAL, El::BigFloat(1),
           d_sdp.free_var_matrix.blocks[block_index], y.blocks[block_index],
           El::BigFloat(-1), dx.blocks[block_index]);
      El::Scale(El::BigFloat(1) / alpha, dx.blocks[block_index]);
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
      El::Scale(El::BigFloat(1) / alpha, dy.blocks[block_index]);
    }
  Block_Vector negative_dx(dx);
  scale_block_vector(negative_dx, El::BigFloat(-1));
  hess_xp.emplace_back(negative_dx, dy);
  // Solve for dx, dy in-place
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal, Q, dx, dy);
  delta_x_y.emplace_back(dx, dy);
}
