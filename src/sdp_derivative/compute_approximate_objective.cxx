#include "../sdp_solve.hxx"

El::BigFloat
compute_approximate_objective(const SDP &sdp, const SDP &d_sdp,
                              const Block_Vector &x, const Block_Vector &y,
                              const Block_Matrix &B_pseudoinverse)
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
      std::cout << "\n " << local_sum;


      // El::Print(d_sdp.primal_objective_c.blocks.at(block),"\ndc");
      // El::Print(d_sdp.free_var_matrix.blocks[block],"dB");
        // El::Print(x.blocks[block],"\nx");
        // El::Print(y.blocks[block],"\ny");
        std::cout << "\n";
        
        // std::cout << "By: "
        //   << By.Height() << " "
        //   << By.Width() << " "
        //   << d_sdp.free_var_matrix.blocks[block].Height() << " "
        //   << d_sdp.free_var_matrix.blocks[block].Width() << " "
        //   << y.blocks[block].Height() << " "
        //   << y.blocks[block].Width() << " "
        //   << "\n";
          
      {
        // dxB = dx.B = db - x.dB
        El::DistMatrix<El::BigFloat> dxB(d_sdp.dual_objective_b);
        El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(-1.0),
                 d_sdp.free_var_matrix.blocks.at(block), x.blocks.at(block),
                 El::BigFloat(1.0), dxB);

        // El::Print(dxB,"\ndxB");
        // std::cout << "\n";
        
        // dxBB_pinv = dx.B.B+
        El::DistMatrix<El::BigFloat> dxBB_pinv(x.blocks.at(block));
        El::Zero(dxBB_pinv);
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 B_pseudoinverse.blocks.at(block), dxB, El::BigFloat(0.0),
                 dxBB_pinv);

        // B_pinv_c = B+.c
        El::DistMatrix<El::BigFloat> B_pinv_c(y.blocks.at(block));
        El::Zero(B_pinv_c);
        El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(1.0),
                 B_pseudoinverse.blocks.at(block),
                 sdp.primal_objective_c.blocks.at(block), El::BigFloat(0.0),
                 B_pinv_c);
        
        // dcBBc = dc - dB.B+.c = dc - db.B_pinv_c
        El::DistMatrix<El::BigFloat> dcBBc(
          d_sdp.primal_objective_c.blocks.at(block));
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(-1.0),
                 d_sdp.free_var_matrix.blocks.at(block),
                 B_pinv_c, El::BigFloat(1.0),
                 dcBBc);


        El::Print(B_pseudoinverse.blocks.at(block),"\nB_pinv");
        El::Print(B_pinv_c,"\nB_pinv_c");
        El::Print(sdp.primal_objective_c.blocks.at(block),"\nc");
        El::Print(d_sdp.primal_objective_c.blocks.at(block),"\ndc");
        
        El::Print(dxBB_pinv,"\ndxBB_pinv");
        El::Print(dcBBc,"\ndcBBc");
        std::cout << "\n";
        
        // Second order term
        // 2*(db - x.dB).B+.(dc - dB.B+.c)
        local_sum += 2*El::Dotu(dxBB_pinv, dcBBc);
      }
    }
  if(!x.blocks.empty() && x.blocks.at(0).Grid().Rank() != 0)
    {
      local_sum = 0;
    }
  return objective
         + El::mpi::AllReduce(local_sum, El::mpi::SUM, El::mpi::COMM_WORLD);
}
