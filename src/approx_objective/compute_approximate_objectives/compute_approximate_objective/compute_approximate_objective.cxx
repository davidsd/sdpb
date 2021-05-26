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

  El::BigFloat linear(0), quad(0);
  // El::BigFloat first(0), second(0);
  // db.y
  objective += El::Dot(d_sdp.dual_objective_b, y.blocks.at(0));

  linear+=El::Dot(d_sdp.dual_objective_b, y.blocks.at(0));
  
  // db.dy/2
  objective += El::Dot(d_sdp.dual_objective_b, dy.blocks.at(0))/2;

  quad += El::Dot(d_sdp.dual_objective_b, dy.blocks.at(0))/2;
  
  El::BigFloat local_sum(0);
  for(size_t block(0); block != x.blocks.size(); ++block)
    {
      // dc.x
      local_sum
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), x.blocks.at(block));

      linear+=Dotu(d_sdp.primal_objective_c.blocks.at(block), x.blocks.at(block));
      
      // dc.dx/2
      local_sum
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), dx.blocks.at(block))/2;

      quad += Dotu(d_sdp.primal_objective_c.blocks.at(block), dx.blocks.at(block))/2;
      {
        // temp = dB.y
        El::DistMatrix<El::BigFloat> temp(x.blocks.at(block));
        El::Zero(temp);
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], y.blocks.at(0), El::BigFloat(0.0),
                 temp);

        // -x.dB.y/
        local_sum -= El::Dotu(temp, x.blocks.at(block));

        linear -= El::Dotu(temp, x.blocks.at(block));
        
        // -dx.dB.y/2
        local_sum -= El::Dotu(temp, dx.blocks.at(block))/2;

        quad -= El::Dotu(temp, dx.blocks.at(block))/2;
        
        // temp = dB.dy
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], dy.blocks.at(0),
                 El::BigFloat(0.0), temp);

        // -x.dB.dy/2
        local_sum -= El::Dotu(temp, x.blocks.at(block))/2;
        
        quad -= El::Dotu(temp, x.blocks.at(block))/2;
      }
    }
  if(!x.blocks.empty() && x.blocks.at(0).Grid().Rank() != 0)
    {
      local_sum = 0;
    }

  return objective
         + El::mpi::AllReduce(local_sum, El::mpi::SUM, El::mpi::COMM_WORLD);
}
