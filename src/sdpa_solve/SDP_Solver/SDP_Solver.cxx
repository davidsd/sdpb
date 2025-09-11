#include "../SDP_Solver.hxx"

namespace Sdpb::Sdpa
{
  SDP_Solver::SDP_Solver(const Solver_Parameters &parameters,
                         const Verbosity &verbosity,
                         const bool &require_initial_checkpoint,
                         const Block_Info &block_info, const El::Grid &grid)
      : x(block_info.primal_dimension, 1),
        X(block_info.block_dimensions, block_info.block_indices, grid),
        Y(X),
        primal_residues(X),
        dual_residues(x),
        current_generation(0)
  {
    if(!load_checkpoint(parameters.checkpoint_in, block_info, verbosity,
                        require_initial_checkpoint))
      {
        X.set_zero();
        Y.set_zero();
        El::Zero(x);

        // X = \Omega_p I
        X.add_diagonal(parameters.initial_matrix_scale_primal);
        // Y = \Omega_d I
        Y.add_diagonal(parameters.initial_matrix_scale_dual);
      }
    if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
      {
        El::Output(
          "Initialize SDP solver"
          "\n\tprimal dimension: ",
          block_info.primal_dimension, // "\n\tdual dimension: ", dual_dimension,
          "\n\tSDP blocks: ", block_info.num_blocks());
      }
  }
}
