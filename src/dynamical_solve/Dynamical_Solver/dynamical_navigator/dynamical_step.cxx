#include <vector>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include "../../../sdp_solve.hxx"
#include "../../../dynamical_solve.hxx"
#include "../../../sdp_read.hxx"
#include "imports.hxx"
#include "../../../approx_objective/Approx_Objective.hxx"

//Compute dx and dy of the central sdp as the standard sdp_solver does in predictor phase. 
//Correspond to - H^-1_xx Del_p L_mu in Eq(13).
//Return: void.
//Update dx, dy. 
void internal_predictor_direction(
  const Block_Info &block_info, const SDP &sdp, const Dynamical_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &grad_x, Block_Vector &grad_y, Block_Vector &dx, Block_Vector &dy, Block_Diagonal_Matrix &R);

//Compute the linear difference of the lagrangians of two sdps 
//The same as the calculation done by Approx_Objective
El::BigFloat compute_delta_lag(const SDP &d_sdp, const Dynamical_Solver &solver);


//Given delta_p(sdp) , compute the (delta_x, delta_y) = H_xx^-1 H_xp delta_p
//as shown in Eq(15). 
//Return: void. 
//Update: hess_xp = H_xp = (RHS(p1)/alpha, RHS(p2)/alpha, ... ), stored to compute the second term on the LHS of Eq(13) 
//        delta_x_y = - H^-1_xx H_xp = - H^-1_xx hess_xp 
void mixed_hess(
  const Block_Info &block_info, const SDP &d_sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q,
  const El::BigFloat &alpha,
  std::vector<std::pair<Block_Vector, Block_Vector>> &hess_xp,
  std::vector<std::pair<Block_Vector, Block_Vector>> &delta_x_y);

  
//Given objectives on a grid in the external parameter (p) space,
//where ePlus correspond to (+e_i), eMinus corresponds to (-e_i) and eSum corresponds to (e_i+e_j) directions respectively, 
//Compute H_pp and Del_p(L_mu) to solve Eq(13).
//Return: void.
//Update: hess = H_pp ,
//        grad = Del_p(L_mu).
void external_grad_hessian(const El::Matrix<El::BigFloat> &ePlus, 
                           const El::Matrix<El::BigFloat> &eMinus, 
                           const El::Matrix<El::BigFloat> &eSum,
                           const El::BigFloat &alpha,
                           El::Matrix<El::BigFloat> &grad,
                           El::Matrix<El::BigFloat> &hess);

//This function is the same as in the usual SDP_Solver, 
//modified to take the argument 'solver' of the right type. 
//Currently not used in the dynamical_step 
void compute_search_direction(
  const Block_Info &block_info, const SDP &sdp, const Dynamical_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const bool &is_corrector_phase, const El::DistMatrix<El::BigFloat> &Q,
  Block_Vector &dx, Block_Diagonal_Matrix &dX, Block_Vector &dy,
  Block_Diagonal_Matrix &dY);

void compute_find_zeros(
  const El::BigFloat &duality_gap, const El::BigFloat &primal_objective,
  bool &find_zeros);

// A subroutine called by run_dynamical
// The external parameter step is passed by the argument 'external_step'
// Compute external_step using Eq (13)
// Compute (dx, dy, dX, dY) using Eq(12)
// Scale both external_step and (dx, dy, dX, dY) by the step length 
void Dynamical_Solver::dynamical_step(
  const Dynamical_Solver_Parameters &dynamical_parameters, 
  const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const El::Grid &grid, 
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y_cholesky,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const Block_Vector &primal_residue_p, El::BigFloat &mu, 
  El::BigFloat &beta, El::BigFloat &primal_step_length, 
  El::BigFloat &dual_step_length, bool &terminate_now, Timers &timers, 
  bool &update_sdp, bool &find_zeros, El::Matrix<El::BigFloat> &external_step) 
{
  auto &step_timer(timers.add_and_start("run.step"));


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
  if(mu > dynamical_parameters.solver_parameters.max_complementarity)
    {
      terminate_now = true;
      return;
    }

  // -----
  //Compute the predictor solution for (dx, dX, dy, dY)
  // -----
    auto &predictor_timer(
      timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));
    beta  = predictor_centering_parameter(dynamical_parameters.solver_parameters, is_primal_and_dual_feasible); 

    // Internal_step: compute dx and dy for the central sdp as in compute_search_direction()      
    //                - H^-1_xx Del_x L_mu in Eq (12) and Eq(13)
    //                Notice that the negative sign has been accounted. 
    Block_Vector internal_dx(x), internal_dy(y);
    Block_Vector grad_x(internal_dx), grad_y(internal_dy);
    Block_Diagonal_Matrix R(X);
    internal_predictor_direction(block_info, sdp, *this, schur_complement_cholesky,
                         schur_off_diagonal, X_cholesky, beta,
                         mu, primal_residue_p, Q, grad_x, grad_y, internal_dx, internal_dy, R); 
   
    //Delta_xy = - H^-1_xx H_xp as the second term on the LHS of Eq(13)
    std::vector<std::pair<Block_Vector, Block_Vector>> H_xp, Delta_xy; 
    int n_external_parameters = dynamical_parameters.n_external_parameters;

    // Eq(13):
    // (H_pp - H_px H^-1_xx H_xp)      dp     =  - Del_p L   + H_px (H^-1_xx Del_x L)
    // (H_pp - (H_px H^-1_xx H_xp))    dp     =  - Del_p L   + H_px . (- internal_dx, - internal_dy)  
    // (hess_pp - hess_mixed)          dp     =  - grad_p    + grad_mixed 

    // External_step: compute hess_mixed := H_px H^-1_xx H_xp   , the second term on the LHS of Eq(13)
    //                        hess_pp    := H_pp                , the first term on the LHS of Eq(13)
    //                        grad_p     := Del_p L             , the first term on the RHS of Eq(13)
    //                        grad_mixed := H_px H^-1_xx Del_x L, the second term on the RHS of Eq(13)
    //                if criteria base on quantities, such as duality_gap, are met
    //                which is decided by function "compute_update_sdp()" called by run_dynamical()
    
    /*
    if (!find_zeros){
       compute_find_zeros(duality_gap, primal_objective, find_zeros);
       if (find_zeros) {std::cout<< "start looking for zeros" << '\n';}
    }
    if (find_zeros && primal_objective < 0 && primal_objective > - 0.0001)
       {std::cout<< "primal_obj "<< primal_objective << '\n'; update_sdp = false;}
    */
    if (update_sdp) 
      {
        El::Matrix<El::BigFloat> grad_p(n_external_parameters,1); // Del_p L
        El::Matrix<El::BigFloat> hess_pp(n_external_parameters, n_external_parameters); //H_pp

        El::Matrix<El::BigFloat> eplus(n_external_parameters,1),eminus(n_external_parameters,1), esum(n_external_parameters, n_external_parameters);
        El::Zero(eplus);
        El::Zero(eminus);
        El::Zero(esum);
        if(dynamical_parameters.new_sdp_path.extension() == ".nsv")
          {
            for(auto &filename : read_file_list(dynamical_parameters.new_sdp_path))
               { 
          	//Assume that the filename takes the form "plus_i","minus_i" and "sum_i_j", f
          	//standing for the change in positive e_i, negative e_i and (e_i + e_j) directions respectively 
                  std::string file_name = filename.stem().string();
                  std::vector<std::string> directions;
                  boost::algorithm::split(directions, file_name, boost::is_any_of("_"));
               	  SDP new_sdp(filename, block_info, grid), d_sdp(new_sdp);
                  Axpy(El::BigFloat(-1), sdp, d_sdp);
                  Approx_Objective approx_obj(block_info, sdp, d_sdp, x, y,
                                           schur_complement_cholesky,
                                           schur_off_diagonal, Q);
                  if (directions[0] == "plus")
                    {
                       mixed_hess(block_info, d_sdp, x, y, schur_complement_cholesky, schur_off_diagonal, Q, dynamical_parameters.alpha, H_xp, Delta_xy);
                       eplus(std::stoi(directions[1])) = approx_obj.d_objective + approx_obj.dd_objective; 
                       //compute_delta_lag(d_sdp, *this);
                    }
                  else if (directions[0] == "minus")
                     {
                       eminus(std::stoi(directions[1])) = approx_obj.d_objective + approx_obj.dd_objective; 
                       //compute_delta_lag(d_sdp,*this); 
                     }
                  else if (directions[0] == "sum")
                     {
                       esum(std::stoi(directions[1]),std::stoi(directions[2])) = approx_obj.d_objective + approx_obj.dd_objective; 
                       //compute_delta_lag(d_sdp,*this); 
                     }
                 }
          }
        else 
          {
             throw std::invalid_argument( "A list of perturbed sdp files are required" );
          }
    
        external_grad_hessian(eplus, eminus, esum, dynamical_parameters.alpha, grad_p, hess_pp);   

        El::Matrix<El::BigFloat> hess_mixed(n_external_parameters,n_external_parameters); //H_px H^-1_xx H_xp in Eq(13).
        El::Matrix<El::BigFloat> grad_mixed(n_external_parameters,1); //H_px (-internal_dx_dy) = H_px (H^-1_xx Del_x L_mu) in Eq (13).  
        for (int i=0; i<n_external_parameters; i++)
          { 
            for (int j=0; j<n_external_parameters; j++)
              {
                // The minus sign compensate the minus sign when calculating Delta_xy in Eq(15)
                hess_mixed(i,j) = - (dot (H_xp.at(i).first,Delta_xy.at(j).first) + dot(H_xp.at(i).second,Delta_xy.at(j).second)); 
              }
                // The minus sign compensates the minus sign when calculating the internal step 
                grad_mixed(i) = - (dot(H_xp.at(i).first,internal_dx) + dot(H_xp.at(i).second,internal_dy));
          }  

        /*
        std::cout<<'\n'
         <<  "eplus: "   << eplus(0,0)<<'\n'
         <<  "eminus: "   << eminus(0,0)<<'\n'
         << "hess_pp " << hess_pp(0,0) <<'\n'
         << "hess_mixed " << hess_mixed(0,0) <<'\n'
         << "Del_p L " << grad_p(0,0) <<'\n'
         << "internal Del L " << grad_mixed(0,0) <<'\n' << std::flush;
        */
 
        // Eq(13):
        // (H_pp - H_px H^-1_xx H_xp)      dp     =  - Del_p L   + H_px (H^-1_xx Del_x L)
        // (H_pp - (H_px H^-1_xx H_xp))    dp     =  - Del_p L   + H_px . (- internal_dx, - internal_dy)  
        // (hess_pp - hess_mixed)          dp     =  - grad_p    + grad_mixed 

        // H_pp - H_px H^-1_xx H_xp, LHS of Eq(13)
        hess_pp -= hess_mixed; 
        // H_px (H^-1_xx Del_x L_mu) - Del_p L_mu, RHS of Eq(13) 
        grad_mixed -= grad_p; 

        external_step = grad_mixed; // we can get rid of grad_mixed here but it might be confusing. 
        El::LinearSolve(hess_pp, external_step);

        /*
        if (find_zeros)
          { 
            std::cout << '\n'
               << "Looking for zeros" << '\n';
            std::cout<< "primal_obj "<< primal_objective << '\n';
            El::Zero(external_step);
            El::BigFloat f, fprime;
            f = dot(grad_x, internal_dx) + dot(grad_y, internal_dy) - primal_objective; 
            fprime = grad_p(0) - (dot(grad_x,Delta_xy.at(0).first) + dot(grad_y,Delta_xy.at(0).second));
            external_step(0) = - f/fprime;
            std::cout<< "step (0): " << external_step(0) << "step size: " << El::Nrm2(external_step) << '\n'; 
            std::cout<< dynamical_parameters.update_sdp_threshold_max << '\n';
          }
        */ 
 
        external_step_size = El::Nrm2(external_step);
        if (external_step_size > dynamical_parameters.update_sdp_threshold_max || external_step_size < dynamical_parameters.update_sdp_threshold_min) 
          { update_sdp = false; } 
     }  


    Block_Vector dx(internal_dx), dy(internal_dy);
    // Update dx and dy if the external parameter step is small enough 
    // RHS of Eq(12)
    //    - H^-1_xx Del_x L_mu - H^-1_xx H_xp dp 
    // =  internal_dx, internal_dy + Delta_xy . external_step
    if (update_sdp) 
      {
        for (int i=0; i<n_external_parameters; i++)
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
      } 
            
    Block_Diagonal_Matrix dX(X), dY(Y);
    // dX = PrimalResidues + \sum_p A_p dx[p]
    constraint_matrix_weighted_sum(block_info, sdp, dx, dX);
    dX += primal_residues;
    
    // dY = Symmetrize(X^{-1} (R - dX Y))
    multiply(dX, Y, dY);
    dY -= R;
    cholesky_solve(X_cholesky, dY);
    dY.symmetrize();
    dY *= El::BigFloat(-1);

    predictor_timer.stop();


  // -----
  // Compute the corrector solution for (dx, dX, dy, dY)
  // -----
    if (external_step_size < dynamical_parameters.update_sdp_threshold_min && external_step_size > 0)
      {
        auto &corrector_timer(
        timers.add_and_start("run.step.computeSearchDirection(betaCorrector)"));
        beta = corrector_centering_parameter(
               dynamical_parameters.solver_parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
               total_psd_rows);
     
        compute_search_direction(block_info, sdp, *this, schur_complement_cholesky,
                                 schur_off_diagonal, X_cholesky,beta,
                                 mu, primal_residue_p, true, Q, dx, dX, dy, dY);
        corrector_timer.stop();
      }


  // Compute step-lengths that preserve positive definiteness of X, Y
  primal_step_length
    = step_length(X_cholesky, dX, dynamical_parameters.solver_parameters.step_length_reduction,
                  "run.step.stepLength(XCholesky)", timers);

  dual_step_length
    = step_length(Y_cholesky, dY, dynamical_parameters.solver_parameters.step_length_reduction,
                  "run.step.stepLength(YCholesky)", timers);

  // If our problem is both dual-feasible and primal-feasible,
  // ensure we're following the true Newton direction.
  if(is_primal_and_dual_feasible)
    {
      primal_step_length = El::Min(primal_step_length, dual_step_length);
      dual_step_length = primal_step_length;
    }
  // step_length = miN
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

  // Update the external step external_step = dualStepLength * external_step
  external_step *= dual_step_length; 
  external_step_size = El::Nrm2(external_step);
  if (external_step_size > dynamical_parameters.update_sdp_threshold_max || external_step_size < dynamical_parameters.update_sdp_threshold_min)
  { update_sdp = false; }

  step_timer.stop();
} 
