#include "../../../sdp_solve.hxx"
#include "../../../dynamical_solve.hxx"
#include <El.hpp>

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

// SDP_solver 
// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
// Lagrangian = new_dual_residues.x + new_b.y + Tr(X,Y) - mu log det (X) 

El::BigFloat compute_delta_lag(const SDP &d_sdp, const Dynamical_Solver &solver){
  El::BigFloat d_objective;
  
  d_objective = 0; 
  
  // dconst
  d_objective += d_sdp.objective_const;
  // db.y 
  d_objective += El::Dot(d_sdp.dual_objective_b, solver.y.blocks.at(0));
  
  El::BigFloat local_diff(0);
  for(size_t block(0); block != solver.x.blocks.size(); ++block)
    {
      // dc.x
      local_diff += Dotu(d_sdp.primal_objective_c.blocks.at(block), solver.x.blocks.at(block));
      // temp = dB.y
      El::DistMatrix<El::BigFloat> temp(solver.x.blocks.at(block));
      El::Zero(temp);
      El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
               d_sdp.free_var_matrix.blocks[block], solver.y.blocks.at(0),
               El::BigFloat(0.0), temp);
      // -x.dB.y/
      local_diff -= El::Dotu(temp, solver.x.blocks.at(block));
   }

    d_objective
    += El::mpi::AllReduce(local_diff, El::mpi::SUM, El::mpi::COMM_WORLD);
  

  return d_objective; 
}

