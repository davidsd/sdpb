#include "../SDP_Solver.hxx"

SDP_Solver::SDP_Solver(const Solver_Parameters &parameters,
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
  if(!load_checkpoint(parameters.checkpoint_in, block_info, verbosity,
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
      X.add_diagonal(parameters.initial_matrix_scale_primal);
      // Y = \Omega_d I
      Y.add_diagonal(parameters.initial_matrix_scale_dual);
    }
  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      const size_t num_blocks = block_info.num_points.size();
      size_t primal_dimension = 0;
      for(size_t i = 0; i < num_blocks; ++i)
        primal_dimension += block_info.get_schur_block_size(i);
      const size_t dual_dimension = dual_objective_b_height;
      El::Output("Initialize SDP solver"
                 "\n\tprimal dimension: ",
                 primal_dimension, "\n\tdual dimension: ", dual_dimension,
                 "\n\tSDP blocks: ", num_blocks);
    }
}
