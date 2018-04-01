#include "../SDP_Solver.hxx"

// Create and initialize an SDPSolver for the given SDP and
// SDP_Solver_Parameters
SDP_Solver::SDP_Solver(const SDP &sdp, const SDP_Solver_Parameters &parameters)
    : sdp(sdp), parameters(parameters), x(sdp.primal_objective.size(), 0),
      X(sdp.psd_matrix_block_dims()), y(sdp.dual_objective.size(), 0), Y(X),
      dx(x), dX(X), dy(y), dY(Y), primal_residues(X), dual_residues(x),
      X_cholesky(X), Y_cholesky(X), Z(X), R(X),
      bilinear_pairings_X_Inv(sdp.bilinear_pairing_block_dims()),
      bilinear_pairings_Y(bilinear_pairings_X_Inv),
      schur_complement(sdp.schur_block_dims()),
      schur_complement_cholesky(schur_complement),
      schur_off_diagonal(sdp.free_var_matrix),
      Q(sdp.free_var_matrix.cols, sdp.free_var_matrix.cols),
      Q_pivots(sdp.free_var_matrix.cols), dy_extended(Q.rows),
      step_matrix_workspace(X)
{
  // initialize bilinearPairingsWorkspace, eigenvaluesWorkspace, QRWorkspace
  for(unsigned int b = 0; b < sdp.bilinear_bases.size(); b++)
    {
      bilinear_pairings_workspace.push_back(
        Matrix(X.blocks[b].rows, bilinear_pairings_X_Inv.blocks[b].cols));
      eigenvalues_workspace.push_back(Vector(X.blocks[b].rows));
      QR_workspace.push_back(Vector(3 * X.blocks[b].rows - 1));
    }

  // X = \Omega_p I
  X.add_diagonal(parameters.initial_matrix_scale_primal);
  // Y = \Omega_d I
  Y.add_diagonal(parameters.initial_matrix_scale_dual);
}
