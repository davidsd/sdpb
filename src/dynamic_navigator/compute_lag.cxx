#include "../sdp_solve.hxx"
#include <El.hpp>

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

// SDP_solver 
// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
// Lagrangian = new_dual_residues.x + new_b.y + Tr(X,Y) - mu log det (X) 

El::BigFloat compute_lag(const SDP &sdp,const SDP &d_sdp, const SDP_Solver &solver){
  El::BigFloat objective;
  
  //b.y
  objective = El::Dot(sdp.dual_objective_b, solver.y.blocks.at(0)) + sdp.objective_const;
  
  // dconst
  objective += d_sdp.objective_const;
  // new_b.y = db.y + b.y 
  objective += El::Dot(d_sdp.dual_objective_b, solver.y.blocks.at(0));
  
  //new_dual_residues.x = new_c.x - x.new_B.y - x.Tr(A_* Y) (A stays the same)
  //  = dual_residues.x + dc.x - x.dB.y 
  for(size_t block(0); block != solver.x.blocks.size(); ++block)
    {
      //dual_residues.x 
      objective += El::Dotu(solver.dual_residues.blocks.at(block), solver.x.blocks.at(block));
      // dc.x
      objective += Dotu(d_sdp.primal_objective_c.blocks.at(block), solver.x.blocks.at(block));
      // temp = dB.y
      El::DistMatrix<El::BigFloat> temp(solver.x.blocks.at(block));
      El::Zero(temp);
      El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
               d_sdp.free_var_matrix.blocks[block], solver.y.blocks.at(0),
               El::BigFloat(0.0), temp);
      // -x.dB.y/
      objective -= El::Dotu(temp, solver.x.blocks.at(block));
   }
  
  //Tr(XY) 
  objective += frobenius_product_symmetric(solver.X, solver.Y);

  //Lagrangian = b.y + dual_residues.x + Tr(XY) - mu log det (X) ?
  //TODO? 
  return objective; 
}

