#include "../sdp_solve.hxx"

void compute_dx_dy(const SDP &d_sdp,
                   const Block_Vector &x, const Block_Vector &y,
                   const Block_Diagonal_Matrix &Y,
                   Block_Vector &dx, Block_Vector &dy);

El::BigFloat
compute_approximate_objective(const SDP &sdp, const SDP &d_sdp,
                              const Block_Vector &x, const Block_Vector &y,
                              const Block_Diagonal_Matrix &Y)
{
  // Constant + first order term: b.y + db.y
  El::BigFloat objective(El::Dot(sdp.dual_objective_b, y.blocks.at(0))
                         + El::Dot(d_sdp.dual_objective_b, y.blocks.at(0)));

  std::cout << "objective: "
            << El::Dot(sdp.dual_objective_b, y.blocks.at(0)) << "\n "
            << El::Dot(d_sdp.dual_objective_b, y.blocks.at(0));

  Block_Vector By(y);
  El::BigFloat local_sum(0);
  for(size_t block(0); block != By.blocks.size(); ++block)
    {
      // First order term: dc.x
      local_sum
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), x.blocks.at(block));

      std::cout << "\n " << local_sum;
      {
        // dBy = dB.y
        El::DistMatrix<El::BigFloat> dBy(x.blocks.at(block));
        El::Zero(dBy);
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], y.blocks[block],
                 El::BigFloat(0.0), dBy);
        // First order term: -x.dB.y
        local_sum -= El::Dotu(dBy, x.blocks[block]);
      }
      std::cout << "\n " << local_sum << "\n";;
    }
  if(!x.blocks.empty() && x.blocks.at(0).Grid().Rank() != 0)
    {
      local_sum = 0;
    }

  Block_Vector dx(x), dy(y);
  // compute_dx_dy(d_sdp, x, y, Y, dx, dy);
  
  return objective
         + El::mpi::AllReduce(local_sum, El::mpi::SUM, El::mpi::COMM_WORLD);
}
