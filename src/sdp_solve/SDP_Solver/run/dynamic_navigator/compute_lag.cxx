#include "../../../../sdp_solve.hxx"
#include <El.hpp>

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

// SDP_solver 
// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
// Lagrangian = new_dual_residues.x + new_b.y + Tr(X,Y) - mu log det (X) 

El::BigFloat compute_lag(const SDP &sdp,const SDP &d_sdp, const SDP_Solver &solver){
  El::BigFloat d_objective;
  
  //b.y
  d_objective = 0; //El::Dot(sdp.dual_objective_b, solver.y.blocks.at(0)) + sdp.objective_const;
  
  // dconst
  d_objective += d_sdp.objective_const;
  // new_b.y = db.y + b.y 
  d_objective += El::Dot(d_sdp.dual_objective_b, solver.y.blocks.at(0));
  std::cout << "dconst + db.y" << ": " << d_objective << '\n' ; 
  //new_dual_residues.x = new_c.x - x.new_B.y - x.Tr(A_* Y) (A stays the same)
  //  = dual_residues.x + dc.x - x.dB.y 
  El::BigFloat local_diff(0);
  for(size_t block(0); block != solver.x.blocks.size(); ++block)
    {
      //dual_residues.x 
      //objective += El::Dotu(solver.dual_residues.blocks.at(block), solver.x.blocks.at(block));
      // dc.x
      local_diff += Dotu(d_sdp.primal_objective_c.blocks.at(block), solver.x.blocks.at(block));
      std::cout << "block" << block << ": " << "dc.x" << local_diff<< '\n' ;
      // temp = dB.y
      El::DistMatrix<El::BigFloat> temp(solver.x.blocks.at(block));
      El::Zero(temp);
      El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
               d_sdp.free_var_matrix.blocks[block], solver.y.blocks.at(0),
               El::BigFloat(0.0), temp);
      // -x.dB.y/
      local_diff -= El::Dotu(temp, solver.x.blocks.at(block));
      std::cout << "block" << block << ": " << "-x.dB.y" << local_diff << '\n' ;
   }

    d_objective
    += El::mpi::AllReduce(local_diff, El::mpi::SUM, El::mpi::COMM_WORLD);
  
  //Tr(XY) 
  //objective += frobenius_product_symmetric(solver.X, solver.Y);

  //Lagrangian = b.y + dual_residues.x + Tr(XY) - mu log det (X) ?
  //TODO? 
  return d_objective; 
}

