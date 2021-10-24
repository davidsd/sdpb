#include <El.hpp>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <Eigen/Dense>
#include "cblas.h"
#include "../sdp_solve.hxx"
#include "../sdp_read.hxx"

typedef Eigen::Matrix<El::BigFloat, Eigen::Dynamic, Eigen::Dynamic> MatrixXBF;
typedef Eigen::Matrix<El::BigFloat, Eigen::Dynamic, 1> VectorXBF;

void internal_predictor_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);

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
  for (unsigned i=0; i<ePlus.Height(); i++)
    { 
      //grad(i) -= eMinus(i);
      //grad(i) /= (2*alpha);
      hess(i,i) += (ePlus(i) + eMinus(i) )/(alpha * alpha);
      for (unsigned j=0; j<ePlus.Height() && i != j ; j++)
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

void setup_solver(const Block_Info &block_info, const El::Grid &grid,
                  const SDP &sdp, const boost::filesystem::path &solution_dir,
                  Block_Diagonal_Matrix &schur_complement_cholesky,
                  Block_Matrix &schur_off_diagonal,
                  El::DistMatrix<El::BigFloat> &Q);
void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);




//We will integrate this into SDP_Solver at some point, 
//then we will not pass in const SDP_Solver &solver.

void dynamic_step(
  const Solver_Parameters &parameters, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const SDP_Solver &solver, const El::Grid &grid, 
  const boost::filesystem::path &input_path, const boost::filesystem::path &solution_dir, 
  //const Block_Diagonal_Matrix &schur_complement_cholesky,
  //const Block_Matrix &schur_off_diagonal, const El::DistMatrix<El::BigFloat> &Q,
  const int &external_dim, const El::BigFloat &alpha,  const Block_Vector &primal_residue_p, El::BigFloat &mu,Timers &timers) 
  
{
     // read x, y, X, Y from file 
      Block_Vector x(block_info.schur_block_sizes(), 
                     block_info.block_indices, block_info.num_points.size(), grid),
                   y(std::vector<size_t>(block_info.num_points.size(), sdp.dual_objective_b.Height()),
                     block_info.block_indices, block_info.num_points.size(), grid);
      Block_Diagonal_Matrix X(block_info.psd_matrix_block_sizes(), 
                              block_info.block_indices, block_info.num_points.size(), grid),
                            Y (X);
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          read_text_block(x.blocks.at(block), solution_dir, "x_",
                          block_index);
          read_text_block(y.blocks.at(block),
                          solution_dir / "y.txt");
          for(size_t psd_block(0); psd_block < 2; ++psd_block)
            {
              if(X.blocks.at(2 * block + psd_block).Height() != 0)
                {
                  const size_t psd_index(2 * block_index + psd_block);
                  read_text_block(X.blocks.at(2 * block + psd_block),
                                  solution_dir, "X_matrix_", psd_index);
                  read_text_block(Y.blocks.at(2 * block + psd_block),
                                  solution_dir, "Y_matrix_", psd_index);
                }
            }
         }

//      //Or get from solver 
//      Block_Vector x(solver.x), y(solver.y); 
//      Block_Diagonal_Matrix X(solver.X), Y(solver.Y);
//
      Block_Diagonal_Matrix X_cholesky(X);
      cholesky_decomposition(X, X_cholesky);

      Block_Diagonal_Matrix schur_complement_cholesky(
        block_info.schur_block_sizes(), block_info.block_indices,
        block_info.num_points.size(), grid);

      Block_Matrix schur_off_diagonal(sdp.free_var_matrix);

      El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
                                     sdp.dual_objective_b.Height());

      setup_solver(block_info, grid, sdp, solution_dir,
                   schur_complement_cholesky, schur_off_diagonal, Q);

      
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

      internal_predictor_direction(block_info, sdp, solver, schur_complement_cholesky,
                           schur_off_diagonal, X_cholesky, beta_predictor,
                           mu, primal_residue_p, Q, internal_dx, internal_dy); 
      
      // approx_step and external Hessian 
      std::vector<std::pair<Block_Vector, Block_Vector>> H_px, Delta_xy;  
//      typedef Eigen::Matrix<El::BigFloat, Eigen::Dynamic, Eigen::Dynamic> MatrixXBF;
//      typedef Eigen::Matrix<El::BigFloat, Eigen::Dynamic, 1> VectorXBF;   
      //std::vector<El::BigFloat> eplus(external_dim), eminus(external_dim);
      //std::vector<std::vector<El::BigFloat>> esum(external_dim, std::vector<El::BigFloat>(external_dim));
      El::Matrix<El::BigFloat> eplus(external_dim,1),eminus(external_dim,1), esum(external_dim, external_dim);
      if(input_path.extension() == ".nsv")
        {
          for(auto &filename : read_file_list(input_path))
             {  
		//Assume that the filename takes the form "plus_i","minus_i" and "sum_i_j", f
		//standing for the change in positive e_i, negative e_i and (e_i + e_j) directions respectively 
                std::string file_name = filename.string();
                std::vector<std::string> directions;
                boost::algorithm::split(directions, file_name, boost::is_any_of("-"));
		//char *word = strtok ( filename.string().c_str(),"-");
		SDP new_sdp(input_path, block_info, grid), d_sdp(new_sdp);
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
     
//         typedef Eigen::Matrix<El::BigFloat, Eigen::Dynamic, Eigen::Dynamic> MatrixXBF;
//         typedef Eigen::Matrix<El::BigFloat, Eigen::Dynamic, 1> VectorXBF;
	El::Matrix<El::BigFloat> approx_hess(external_dim,external_dim);
        El::Matrix<El::BigFloat> internal_grad(external_dim,1);
//         MatrixXBF approx_hess(external_dim,external_dim);
//         VectorXBF internal_grad(external_dim);
 
         El::Matrix<El::BigFloat>  grad(eplus), hess(esum);
//       std::vector<El::BigFloat> grad(eplus)//, internal_grad(eplus);
//       std::vector<std::vector<El::BigFloat>> hess(esum)//, approx_hess(esum);
//       El::BigFloat *internal_grad, *approx_hess;
//       approx_hess = (El::BigFloat *) malloc(sizeof(El::BigFloat)*external_dim*external_dim);
//       internal_grad = (El::BigFloat *) malloc(sizeof(El::BigFloat)*external_dim); 
       external_grad_hessian(eplus, eminus, esum, alpha, grad, hess);   
        
       for (unsigned i=0; i<esum.Height(); i++)
         { 
           for (unsigned j=0; j<esum.Width(); j++)
             {
               approx_hess(i,j) = dot (H_px.at(i).first,Delta_xy.at(j).first) + dot(H_px.at(i).second,Delta_xy.at(j).second); 
//               approx_hess[external_dim * i+j] = dot (H_px.at(i).first,Delta_xy.at(j).first) + dot(H_px.at(i).second,Delta_xy.at(j).second); 
             }
               internal_grad(i) = dot(H_px.at(i).first,internal_dx) + dot(H_px.at(i).second,internal_dy);
//           internal_grad[i] = dot(H_px.at(i).first,internal_dx) + dot(H_px.at(i).second,internal_dy);
         }  
//      double *a;
//      a = (double *)malloc(sizeof(double)*external_dim*external_dim);
//      for (unsigned i=0; i<external_dim; i++)
//	{
//	for (unsigned j=0; j<external_dim; j++)
//	  {
//	    a[external_dim * j + i] = external_dim * j + i;
//	  }
//        }
//         El::LinearSolve( A, b );	
//       hess - approx_hess 
//       grad - internal_grad
//       El::LinearSolve(hess - approx_hess, grad - internal_grad)
         hess -= approx_hess;
         grad -= internal_grad;
         El::Matrix<El::BigFloat> step(grad);
         El::LinearSolve(hess,step);
} 

