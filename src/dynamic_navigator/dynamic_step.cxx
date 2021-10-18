#include "Approx_Objective.hxx"
#include "../sdp_solve.hxx"



void internal_predictor_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);




//Give delta_sdp, compute (delta_x, delta_y) = H_xx^-1 H_xp delta_p
//Return delta_c_b_B =  H_xp delta_p , delta_x_y = H_xx^-1 H_xp delta_p
approx_step(
  const Block_Info &block_info, const SDP &d_sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q
  std:vector<Block_Vector> &delta_c_b_B,
  std:vector<Block_Vector> &delta_x_y)
{
  Block_Vector dx(x), dy(y);
  
  for(size_t block_index(0); block_index != dx.blocks.size(); ++block_index)
    {
      // dx = -dc + dB.y
      dx.blocks[block_index] = d_sdp.primal_objective_c.blocks[block_index];
      Gemv(El::Orientation::NORMAL, El::BigFloat(1),
           d_sdp.free_var_matrix.blocks[block_index], y.blocks[block_index],
           El::BigFloat(-1), dx.blocks[block_index]);

      // dy = db - x.dB
      El::Zero(dy.blocks[block_index]);
      if(block_info.block_indices[block_index] == 0)
        {
          dy.blocks[block_index] = d_sdp.dual_objective_b;
        }
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               d_sdp.free_var_matrix.blocks[block_index],
               x.blocks[block_index], El::BigFloat(1.0),
               dy.blocks[block_index]);
    }
  delta_c_b_B.push_back(dx, dy);  
  // Solve for dx, dy in-place
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal, Q, dx, dy);
  delta_x_y.push_back(dx, dy);

}
  
//Given objectives on a grid in the external parameter (p) space,
//where ePlus correspond to (+e_i), eMinus corresponds to (-e_i) and eSum corresponds to (e_i+e_j) directions respectively 
//Compute H_pp and grad_p(Lag)
external_grad_hessian(std::vector<El::BigFloat> &ePlus, 
                      std::vector<El::BigFloat> &eMinus, 
                      std::vector<std::vector<El::BigFloat>> &eSum
                      El::BigFloat &alpha)
{
  std::vector<El::BigFloat> grad(ePlus);
  std::vector<std::vector<El::BigFloat>> hess(eSum);
  grad -=  eMinus;
  grad /= (2*alpha);
  for (unsigned i=0; i<ePlus.size(); i++)
    { 
      hess.at(i).at(i) += (ePlus.at(i) + eMinus.at(i) )/pow(alpha,2);
      for (unsigned j=0; j<ePlus.size() && i != j ; j++)
        {
          hess.at(i).at(j) += (- ePlus.at(i) - ePlus.at(j))/pow(alpha,2);
        }
    }
}

//We will integrate this into SDP_Solver at some point, 
//then we will not pass in const SDP_Solver &solver.

dynamic_step(
  const Solver_Parameters &parameters, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const SDP_Solver &solver, const El::Grid &grid, 
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal, const El::DistMatrix<El::BigFloat> &Q,
  const Block_Vector &primal_residue_p, El::BigFloat &mu,Timers &timers) 

{
  El::Environment env(argc, argv);

  try
    {
      Approx_Parameters parameters(argc, argv);
      if(!parameters.is_valid())
        {
          return 0;
        }
      El::gmp::SetPrecision(parameters.precision);
      Block_Info block_info(parameters.sdp_path, parameters.solution_dir,
                            parameters.procs_per_node,
                            parameters.proc_granularity, Verbosity::none);

      El::Grid grid(block_info.mpi_comm.value);
      SDP sdp(parameters.sdp_path, block_info, grid);

      std::vector<size_t> block_offsets(sdp.free_var_matrix.blocks.size() + 1,
                                        0);
      for(size_t p(0); p < sdp.free_var_matrix.blocks.size(); ++p)
        {
          block_offsets[p + 1]
            = block_offsets[p] + sdp.free_var_matrix.blocks[p].Height();
        }

      

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
          read_text_block(x.blocks.at(block), parameters.solution_dir, "x_",
                          block_index);
          read_text_block(y.blocks.at(block),
                          parameters.solution_dir / "y.txt");
          for(size_t psd_block(0); psd_block < 2; ++psd_block)
            {
              if(X.blocks.at(2 * block + psd_block).Height() != 0)
                {
                  const size_t psd_index(2 * block_index + psd_block);
                  read_text_block(X.blocks.at(2 * block + psd_block),
                                  parameters.solution_dir, "X_matrix_", psd_index);
                  read_text_block(Y.blocks.at(2 * block + psd_block),
                                  parameters.solution_dir, "Y_matrix_", psd_index);
                }
            }
         }

      //Or get from solver 

      auto &frobenius_timer(
        timers.add_and_start("run.step.frobenius_product_symmetric"));
      mu = frobenius_product_symmetric(X, Y) / total_psd_rows;
      frobenius_timer.stop();
      if(mu > parameters.max_complementarity)
        {
          terminate_now = true;
          return;
        }

      auto &predictor_timer(
        timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));
      beta_predictor = predictor_centering_parameter(parameters, is_primal_and_dual_feasible); 

      Block_Diagonal_Matrix X_cholesky(X);
      cholesky_decomposition(X, X_cholesky);

      Block_Diagonal_Matrix schur_complement_cholesky(
        block_info.schur_block_sizes(), block_info.block_indices,
        block_info.num_points.size(), grid);
      Block_Matrix schur_off_diagonal(sdp.free_var_matrix);
      El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
                                     sdp.dual_objective_b.Height());
      
      setup_solver(block_info, grid, sdp, parameters.solution_dir,
                   schur_complement_cholesky, schur_off_diagonal, Q);
      
      Block_Vector dx(x), dy(y)

      internal_predictor_direction(block_info, sdp, *this, schur_complement_cholesky,
                           schur_off_diagonal, X_cholesky, beta_predictor,
                           mu, primal_residue_p, Q, dx, dy); 
      

      std:vector<Block_Vector> H_px, Delta_xy;     
      std::vector<std::vector<El::BigFloat>> esum(dim, vector<El::BigFloat>(dim)); approx_hess(esum); 
      if(input_path.extension() == ".nsv")
        {
          for(auto &filename : read_file_list(input_path))
             {  
		//Assume that the filename takes the form "plus_i","minus_i" and "sum_i_j", 
		//standing for the change in positive e_i, negative e_i and (e_i + e_j) directions respectively 
                char *word = strtok ( filename.string,"-");
		SDP new_sdp(input_path, block_info, grid), d_sdp(new_sdp);
                Axpy(El::BigFloat(-1), sdp, d_sdp);
                if (word = "plus")
                  {
                     approx_step(block_info, d_sdp, x, y, schur_complement_cholesky, schur_off_diagonal, Q, H_px, Delta_xy)
                     eplus.at(int(strtok(NULL, " "))) = Lagrangian  
                  }
                else if (word = "minus")
                   {
                     eminus.at(int(strtok(NULL, " "))) = Approx_Objective(block_info, sdp, d_sdp, x, y,
                                           schur_complement_cholesky,
                                           schur_off_diagonal, Q).objective 
                   }
                else if (word = "sum")
                   {
                   esum.at(int(strtok(NULL, " "))).at(int(strtok(NULL, " "))) = Approx_Objective(block_info, sdp, d_sdp, x, y,
                                           schur_complement_cholesky,
                                           schur_off_diagonal, Q).objective 
	           }
               }
        }
       grad,hess = external_grad_hessian(eplus, eminus, esum, alpha)   
        
       for (unsigned i=0; i<eplus.size(); i++)
         { 
           for (unsigned j=0; j<eplus.size(); i++)
             {
               approx_hess.at(i).at(j) = El::Dotu (H_px.at(i),M_delta_x_y.at(j)); 
             }
           internal_grad.at(i) = El::Dotu (H_px.at(i),internal_step)
         }  
       hess - approx_hess 
       grad - internal_grad
       El::LinearSolve(hess - approx_hess, grad - internal_grad)
} 

