#include "../SDP_Solver.hxx"

// Create and initialize an SDPSolver for the given SDP and
// SDP_Solver_Parameters
SDP_Solver::SDP_Solver(const std::vector<boost::filesystem::path> &sdp_files,
                       const SDP_Solver_Parameters &parameters)
    : sdp(sdp_files), parameters(parameters),
      x(sdp.primal_objective.size(), 0), X(sdp.psd_matrix_block_dims()),
      y(sdp.dual_objective.size(), 0), Y(X), dx(x), dX(X), dy(y), dY(Y),
      primal_residues(X), dual_residues(x)
{
  // X = \Omega_p I
  X.add_diagonal(parameters.initial_matrix_scale_primal);
  // Y = \Omega_d I
  Y.add_diagonal(parameters.initial_matrix_scale_dual);
}
