#include <El.hpp>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include "cblas.h"
#include "../../../../sdp_solve.hxx"
#include "../../../../sdp_read.hxx"


void internal_predictor_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy, Block_Diagonal_Matrix &R);

El::BigFloat compute_lag(const SDP &sdp,const SDP &d_sdp, const SDP_Solver &solver);


//Give delta_sdp, compute (delta_x, delta_y) = H_xx^-1 H_xp delta_p
//Return delta_c_b_B =  H_xp delta_p , delta_x_y = H_xx^-1 H_xp delta_p
void approx_step(
  const Block_Info &block_info, const SDP &d_sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q,
  std::vector<std::pair<Block_Vector, Block_Vector>> &delta_c_b_B,
  std::vector<std::pair<Block_Vector, Block_Vector>> &delta_x_y);

  
//Given objectives on a grid in the external parameter (p) space,
//where ePlus correspond to (+e_i), eMinus corresponds to (-e_i) and eSum corresponds to (e_i+e_j) directions respectively 
//Compute H_pp and grad_p(Lag)
void external_grad_hessian(const El::Matrix<El::BigFloat> &ePlus, 
                           const El::Matrix<El::BigFloat> &eMinus, 
                           const El::Matrix<El::BigFloat> &eSum,
                           const El::BigFloat &alpha,
                           El::Matrix<El::BigFloat> &grad,
                           El::Matrix<El::BigFloat> &hess)
{
  grad -= eMinus;
  grad *= 1.0/(2.0*alpha);
  for (int i=0; i<ePlus.Height(); i++)
    { 
      hess(i,i) += (ePlus(i) + eMinus(i) )/(alpha * alpha);
      for (int j=0; j<ePlus.Height() && i != j ; j++)
        {
          hess(i,j) += (- ePlus(i) - ePlus(j))/(alpha * alpha);
        }
    }
}


El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

El::BigFloat predictor_centering_parameter(const Solver_Parameters &parameters,
                                           const bool is_primal_dual_feasible);

void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);

void cholesky_solve(const Block_Diagonal_Matrix &ACholesky,
                    Block_Diagonal_Matrix &X);

void constraint_matrix_weighted_sum(const Block_Info &block_info,
                                    const SDP &sdp, const Block_Vector &a,
                                    Block_Diagonal_Matrix &result);

El::BigFloat step_length(const Block_Diagonal_Matrix &MCholesky,
                         const Block_Diagonal_Matrix &dM,
                         const El::BigFloat &gamma,
                         const std::string &timer_name,
                         Timers &timers);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const El::Grid &block_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q,
  Timers &timers);

void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);

// C := alpha*A*B + beta*C
void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C);

// C := A B
inline void multiply(const Block_Diagonal_Matrix &A,
                     const Block_Diagonal_Matrix &B, Block_Diagonal_Matrix &C)
{
  scale_multiply_add(El::BigFloat(1), A, B, El::BigFloat(0), C);
}



//We will integrate this into SDP_Solver at some point, 
//then we will not pass in const SDP_Solver &solver.

void dynamic_step(
  const Solver_Parameters &parameters, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, SDP_Solver &solver, const El::Grid &grid, 
  const boost::filesystem::path &new_sdp_path, 
  const int &n_external_paramters, const El::BigFloat &alpha,  
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y_cholesky,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const Block_Vector &primal_residue_p, 
  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length, El::BigFloat &mu,Timers &timers) 
{
     Block_Vector& x = solver.x, y = solver.y;
      Block_Diagonal_Matrix& X = solver.X, Y = solver.Y;

      Block_Diagonal_Matrix schur_complement_cholesky(
        block_info.schur_block_sizes(), block_info.block_indices,
        block_info.num_points.size(), grid);

      Block_Matrix schur_off_diagonal(sdp.free_var_matrix);

      El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
                                     sdp.dual_objective_b.Height());

      initialize_schur_complement_solver(block_info, sdp, A_X_inv, A_Y, grid,
                                         schur_complement_cholesky,
                                         schur_off_diagonal, Q, timers);
     
      auto &frobenius_timer(
        timers.add_and_start("run.step.frobenius_product_symmetric"));
      mu = frobenius_product_symmetric(X, Y) / total_psd_rows;
      frobenius_timer.stop();
// TODO
//      if(mu > parameters.max_complementarity)
//        {
//          terminate_now = true;
//          return;
//        }

      auto &predictor_timer(
        timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));
      El::BigFloat beta_predictor;
      beta_predictor = predictor_centering_parameter(parameters, is_primal_and_dual_feasible); 

      // Internal_step      
      Block_Vector internal_dx(x), internal_dy(y);
      Block_Diagonal_Matrix R(X);
      internal_predictor_direction(block_info, sdp, solver, schur_complement_cholesky,
                           schur_off_diagonal, X_cholesky, beta_predictor,
                           mu, primal_residue_p, Q, internal_dx, internal_dy, R); 
      
      // approx_step and external Hessian 
      std::vector<std::pair<Block_Vector, Block_Vector>> H_px, Delta_xy;  
      El::Matrix<El::BigFloat> eplus(n_external_paramters,1),eminus(n_external_paramters,1), esum(n_external_paramters, n_external_paramters);
      if(new_sdp_path.extension() == ".nsv")
        {
          for(auto &filename : read_file_list(new_sdp_path))
             {  
		//Assume that the filename takes the form "plus_i","minus_i" and "sum_i_j", f
		//standing for the change in positive e_i, negative e_i and (e_i + e_j) directions respectively 
                std::string file_name = filename.string();
                std::vector<std::string> directions;
                boost::algorithm::split(directions, file_name, boost::is_any_of("-"));
		SDP new_sdp(new_sdp_path, block_info, grid), d_sdp(new_sdp);
                Axpy(El::BigFloat(-1), sdp, d_sdp);
                if (directions[0] == "plus")
                  {
                     approx_step(block_info, d_sdp, x, y, schur_complement_cholesky, schur_off_diagonal, Q, H_px, Delta_xy);
                     eplus(std::stoi(directions[1])) = compute_lag(sdp,d_sdp,solver);
                  }
                else if (directions[0] == "minus")
                   {
                     eminus(std::stoi(directions[1])) = compute_lag(sdp,d_sdp,solver); 
                   }
                else if (directions[0] == "sum")
                   {
                     esum(std::stoi(directions[1]),std::stoi(directions[2])) = compute_lag(sdp,d_sdp,solver); 
	           }
               }
         }
       else 
         {
           throw std::invalid_argument( "A list of perturbed sdp files are required" );
         }
     
 
        El::Matrix<El::BigFloat>  grad(eplus), hess(esum);
        external_grad_hessian(eplus, eminus, esum, alpha, grad, hess);   

        El::Matrix<El::BigFloat> approx_hess(n_external_paramters,n_external_paramters);
        El::Matrix<El::BigFloat> internal_grad(n_external_paramters,1);        
        for (int i=0; i<esum.Height(); i++)
          { 
            for (int j=0; j<esum.Width(); j++)
              {
                approx_hess(i,j) = dot (H_px.at(i).first,Delta_xy.at(j).first) + dot(H_px.at(i).second,Delta_xy.at(j).second); 
              }
                internal_grad(i) = dot(H_px.at(i).first,internal_dx) + dot(H_px.at(i).second,internal_dy);
          }  

         hess -= approx_hess;
         internal_grad -= grad;
         El::Matrix<El::BigFloat> external_step(internal_grad);
         El::LinearSolve(hess,external_step);

         Block_Vector dx(internal_dx), dy(internal_dy);
         for (int i=0; i<esum.Height(); i++)
           {
             for(size_t block = 0; block < x.blocks.size(); ++block)
               {
                 El::Axpy(external_step(i), Delta_xy.at(i).first.blocks[block], dx.blocks[block]);
               }
             for(size_t block = 0; block < dy.blocks.size(); ++block)
               {
                 El::Axpy(external_step(i), Delta_xy.at(i).second.blocks[block], dy.blocks[block]);
               }
          }
 
         
         Block_Diagonal_Matrix dX(X), dY(Y);
         // dX = PrimalResidues + \sum_p A_p dx[p]
         constraint_matrix_weighted_sum(block_info, sdp, dx, dX);
         dX += solver.primal_residues;
         
         // dY = Symmetrize(X^{-1} (R - dX Y))
         multiply(dX, solver.Y, dY);
         dY -= R;
         cholesky_solve(X_cholesky, dY);
         dY.symmetrize();
         dY *= El::BigFloat(-1);
            
         // Compute step-lengths that preserve positive definiteness of X, Y
         primal_step_length
           = step_length(X_cholesky, dX, parameters.step_length_reduction,
                         "run.step.stepLength(XCholesky)", timers);
       
         dual_step_length
           = step_length(Y_cholesky, dY, parameters.step_length_reduction,
                         "run.step.stepLength(YCholesky)", timers);

         // If our problem is both dual-feasible and primal-feasible,
         // ensure we're following the true Newton direction.
         if(is_primal_and_dual_feasible)
           {
             primal_step_length = El::Min(primal_step_length, dual_step_length);
             dual_step_length = primal_step_length;
           }
       
         // Update the primal point (x, X) += primalStepLength*(dx, dX)
         for(size_t block = 0; block < x.blocks.size(); ++block)
           {
             El::Axpy(primal_step_length, dx.blocks[block], x.blocks[block]);
           }
         dX *= primal_step_length;
       
         X += dX;
       
         // Update the dual point (y, Y) += dualStepLength*(dy, dY)
         for(size_t block = 0; block < dy.blocks.size(); ++block)
           {
             El::Axpy(dual_step_length, dy.blocks[block], y.blocks[block]);
           }
         dY *= dual_step_length;
       
         Y += dY;

} 
