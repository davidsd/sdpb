#include "../../../sdp_solve.hxx"
#include "../../../dynamical_solve.hxx"
#include <El.hpp>
#include <math.h>
#include <iostream>
//#include <gmpxx.h>
//#include <boost/multiprecision/gmp.hpp>
//#include <mpfr.h>
//#include <mpf2mpfr.h>
//using namespace boost::multiprecision;

//Does the same thing as Approx_Objective computes the linear difference 

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);


El::Base<El::BigFloat> det_log_cholesky(const Block_Diagonal_Matrix &L){
   El::Base<El::BigFloat> local_det_log(0);
   for(size_t b = 0; b < L.blocks.size(); b++)
     {
       El::SafeProduct<El::Base<El::BigFloat>> safeDet = El::hpd_det::AfterCholesky(El::UpperOrLowerNS::LOWER,L.blocks[b]);
       local_det_log += safeDet.kappa*safeDet.n;
       //Notice that El::hpd_det::AfterCholesky multiply the det by 2 so the result is det(X = LL^T) instead of det(L);
     }
   //std::cout << "before all reduce: " << local_det_log << std::endl;
   El::Base<El::BigFloat> final_result = El::mpi::AllReduce(local_det_log, El::mpi::SUM, El::mpi::COMM_WORLD);
   //std::cout << "after all reduce: " << final_result   << std::endl;

   return final_result;
  
}

// solver 
// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
// Lagrangian = dual_residues.x + b.y + Tr(X,Y) - mu log det (X) 

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

// solver 
// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
// Lagrangian = dual_residues.x + b.y + Tr(X,Y) - mu log det (X) 

El::BigFloat compute_lag(const El::BigFloat mu, const Block_Diagonal_Matrix &X_cholesky, const Dynamical_Solver &solver){
  El::BigFloat lag(0);
  // dual_objective = f + b . y
  lag = solver.dual_objective;
  
  El::BigFloat trXY = frobenius_product_symmetric(solver.X, solver.Y);

  // Tr(XY)
  lag += trXY;

  
  // mu log det X = mu log det LL^T = 2 mu log det L; 
  // The El::HPDDeterminant routine takes log of the diagonal entries and then exponentiate the sum
  // depending on if we want to use the whole routine or just use the sum of the logs  
  //Notice that El::hpd_det::AfterCholesky already accounts for the factor of 2 

  El::BigFloat muLogDetX = mu * (det_log_cholesky(X_cholesky));
  lag -= muLogDetX;

  El::BigFloat local_residues(0);
  for(size_t block(0); block != solver.x.blocks.size(); ++block)
    {
      // dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
      // Lagrangian += dual_residues.x 
      local_residues += Dotu(solver.dual_residues.blocks.at(block), solver.x.blocks.at(block));
   }

  El::BigFloat dual_residue_dot_x = El::mpi::AllReduce(local_residues, El::mpi::SUM, El::mpi::COMM_WORLD);

   lag 
    += dual_residue_dot_x;

   if (El::mpi::Rank() == 0) std::cout << "finite mu nvg :"
	   << " d.x = " << dual_residue_dot_x
	   << " b.y = " << solver.dual_objective
	   << " tr(XY) = " << trXY
	   << " mu = " << mu
	   << " mu*log(detX) = " << muLogDetX
	   << " nvg = " << lag
	   << '\n' << std::flush;

  return lag;
}


//int main(){
//  
//  El::BigFloat lag(1.73, 200);
//  std::cout << "lag: " << lag << " precision " << lag.Precision() << std::endl;
//  mpfr_t float_lag, log_lag;
//  mpfr_set_default_prec(lag.Precision());
//  mpfr_init_set_f(float_lag, lag.gmp_float.get_mpf_t(), MPFR_RNDN);
//  mpfr_init(log_lag);
//  mpfr_log(log_lag, float_lag,MPFR_RNDN);
//  std::cout << "log_lag: " << log_lag << " precision " << mpfr_get_default_prec() << std::endl;
//  El::BigFloat log_BF;
//  //mpf_class mpf_log;
//  mpfr_get_f(log_BF.gmp_float.get_mpf_t(), log_lag, MPFR_RNDN);
//  std::cout << "log BF: " << log_BF << " precision " << log_BF.Precision() << std::endl;
//

  //mpf_float float_lag;
  //mpf_float::default_precision(lag.Precision());

  //std::cout << "lag: " << lag << " precision " << lag.Precision() << std::endl;  
  //mpfr_t float_lag, log_lag;
  //mpfr_set_default_prec(lag.Precision());
  //mpfr_init_set_f(float_lag, lag.gmp_float.get_mpf_t(), MPFR_RNDN);
  //mpfr_init(log_lag);
  //mpfr_log(log_lag, float_lag,MPFR_RNDN);
  //std::cout << "log_lag: " << log_lag << " precision " << mpfr_get_default_prec() << std::endl;
  //El::BigFloat log_BF;
  ////mpf_class mpf_log;
  //mpfr_get_f(log_BF.gmp_float.get_mpf_t(), log_lag, MPFR_RNDN);
  //std::cout << "log BF: " << log_BF << " precision " << log_BF.Precision() << std::endl;
  //El::BigFloat log_BF
  //(mpf_log, lag.Precision());
  //mpf_float 
  //std::cout << "c++ std log" << log(float_lag) << std::endl; 
  // mu log det (X) = 2 mu log det (X_cholesky)
  //El::BigFloat det_X_cholesky(det_block_cholesky(X_cholesky));
  //std::cout << "after det" << std::endl;


//  return 0;
//
//}

