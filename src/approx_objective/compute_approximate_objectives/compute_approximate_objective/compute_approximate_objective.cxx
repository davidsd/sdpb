#include "../../../sdp_solve.hxx"

void compute_dx_dy(const Block_Info &block_info, const SDP &d_sdp,
                   const Block_Vector &x, const Block_Vector &y,
                   const Block_Diagonal_Matrix &schur_complement_cholesky,
                   const Block_Matrix &schur_off_diagonal,
                   const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx,
                   Block_Vector &dy);

El::BigFloat compute_approximate_objective(
  const Block_Info &block_info, const SDP &sdp, const SDP &d_sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q)
{
  Block_Vector dx(x), dy(y);
  compute_dx_dy(block_info, d_sdp, x, y, schur_complement_cholesky,
                schur_off_diagonal, Q, dx, dy);

  // b.y
  El::BigFloat objective(El::Dot(sdp.dual_objective_b, y.blocks.at(0)));
  // ydy = y + dy/2
  El::DistMatrix<El::BigFloat> ydy(y.blocks.at(0));
  El::Axpy(El::BigFloat(0.5), dy.blocks.at(0), ydy);

  // db.(y + dy/2)
  objective += El::Dot(d_sdp.dual_objective_b, ydy);

  El::BigFloat local_sum(0);
  for(size_t block(0); block != x.blocks.size(); ++block)
    {
      // dc.x
      local_sum
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), x.blocks.at(block));
      // dc.dx/2
      local_sum
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), dx.blocks.at(block))
           / 2;
      {
        // temp = dB.(y + dy/2)
        El::DistMatrix<El::BigFloat> temp(x.blocks.at(block));
        El::Zero(temp);
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], ydy, El::BigFloat(0.0),
                 temp);

        // xdx = x + dx/2
        El::DistMatrix<El::BigFloat> xdx(x.blocks[block]);
        El::Axpy(El::BigFloat(0.5), dx.blocks[block], xdx);

        // -(x + dx/2).dB.(y + dy/2)
        // Doing it this way reduces the number of dot products, but
        // includes the incorrect third order term dx.dB.dy/4 which
        // comes along for the ride.
        local_sum -= El::Dotu(temp, xdx);
      }
    }
  if(!x.blocks.empty() && x.blocks.at(0).Grid().Rank() != 0)
    {
      local_sum = 0;
    }

  return objective
         + El::mpi::AllReduce(local_sum, El::mpi::SUM, El::mpi::COMM_WORLD);
}
