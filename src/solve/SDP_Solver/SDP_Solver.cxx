#include "../SDP_Solver.hxx"

// Create and initialize an SDPSolver for the given SDP and
// SDP_Solver_Parameters
SDP_Solver::SDP_Solver(const std::vector<boost::filesystem::path> &sdp_files,
                       const SDP_Solver_Parameters &parameters)
    : sdp(sdp_files), parameters(parameters), x(sdp.schur_block_dims),
      X(sdp.psd_matrix_block_dims()), y(sdp.dual_objective_b.Height(), 1),
      Y(X), primal_residues(X),
      // FIXME: Maybe we can use schur_block_dims() instead of making a copy?
      dual_residues(x)
{
  X.set_zero();
  Y.set_zero();
  for(auto &b : x.blocks)
    {
      Zero(b);
    }
  Zero(y);

  // X = \Omega_p I
  X.add_diagonal(parameters.initial_matrix_scale_primal);
  // Y = \Omega_d I
  Y.add_diagonal(parameters.initial_matrix_scale_dual);
}
