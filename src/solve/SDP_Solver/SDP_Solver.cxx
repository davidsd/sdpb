#include "../SDP_Solver.hxx"

// Create and initialize an SDPSolver for the given SDP and
// SDP_Solver_Parameters
SDP_Solver::SDP_Solver(const std::vector<boost::filesystem::path> &sdp_files,
                       const SDP_Solver_Parameters &parameters)
    : sdp(sdp_files), parameters(parameters),
      x(sdp.primal_objective_c.size(), 0),
      x_elemental(sdp.schur_block_dims()), X(sdp.psd_matrix_block_dims()),
      y(sdp.dual_objective_b_elemental.Height(), 0),
      y_elemental(sdp.dual_objective_b_elemental.Height(), 1), Y(X),
      primal_residues(X), dual_residues(x),
      // FIXME: Maybe we can use schur_block_dims() instead of making a copy?
      dual_residues_elemental(x_elemental)
{
  X.set_zero();
  Y.set_zero();
  for(auto &b : x_elemental.blocks)
    {
      Zero(b);
    }
  Zero(y_elemental);

  // X = \Omega_p I
  X.add_diagonal(parameters.initial_matrix_scale_primal_elemental);
  // Y = \Omega_d I
  Y.add_diagonal(parameters.initial_matrix_scale_dual_elemental);
}
