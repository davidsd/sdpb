#include "../Dynamical_Solver.hxx"

Dynamical_Solver::Dynamical_Solver(const Dynamical_Solver_Parameters &dynamical_parameters,
                       const Verbosity &verbosity,
                       const bool &require_initial_checkpoint,
                       const Block_Info &block_info, const El::Grid &grid,
                       const size_t &dual_objective_b_height)
    : x(block_info.schur_block_sizes(), block_info.block_indices,
        block_info.num_points.size(), grid),
      X(block_info.psd_matrix_block_sizes(), block_info.block_indices,
        block_info.num_points.size(), grid),
      y(std::vector<size_t>(block_info.num_points.size(),
                            dual_objective_b_height),
        block_info.block_indices, block_info.num_points.size(), grid),
      Y(X), primal_residues(X),
      dual_residues(block_info.schur_block_sizes(), block_info.block_indices,
                    block_info.num_points.size(), grid),
      current_generation(0)
{
  if(!load_checkpoint(dynamical_parameters.solver_parameters.checkpoint_in, block_info, verbosity,
                      require_initial_checkpoint))
    {
      X.set_zero();
      Y.set_zero();
      for(auto &block : x.blocks)
        {
          El::Zero(block);
        }
      for(auto &block : y.blocks)
        {
          El::Zero(block);
        }

      // X = \Omega_p I
      X.add_diagonal(dynamical_parameters.solver_parameters.initial_matrix_scale_primal);
      // Y = \Omega_d I
      Y.add_diagonal(dynamical_parameters.solver_parameters.initial_matrix_scale_dual);
    }
  
  
  hess_BFGS_updateQ = false;
  p_step = 1;
  d_step = 1;
  specified_ext_param_Q = false;

}
