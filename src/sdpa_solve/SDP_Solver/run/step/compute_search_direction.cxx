#include "sdpa_solve/SDP_Solver/run/constraint_matrix_weighted_sum.hxx"
#include "sdpa_solve/SDP_Solver.hxx"

// Compute the search direction (dx, dX, dy, dY) for the predictor and corrector phases.

namespace Sdpb::Sdpa
{
  // dx[p] = -dual_residues[i] - Tr(F_i Z)
  void
  compute_schur_RHS(const SDP &sdp, const Primal_Dist_Vector &dual_residues,
                    const Block_Diagonal_Matrix &Z, Primal_Dist_Vector &dx)
  {
    for(size_t p = 0; p < sdp.primal_dimension(); ++p)
      {
        const auto d_p = dual_residues.Get(p, 0);
        const auto tr_F_p_Z = dotu(sdp.sdp_blocks_F.at(p), Z);
        dx.Set(p, 0, -d_p - tr_F_p_Z);
      }
  }

  // Solve:
  // S_ij dx_j = - d_i - Tr(Z F_i)
  // dx passed to this function already contains the RHS (computed in compute_schur_RHS)
  void solve_schur_complement_equation(const El::DistMatrix<El::BigFloat> &S,
                                       Primal_Dist_Vector &dx)
  {
    El::DistMatrix<El::BigFloat> dx_dist = dx;
    // dx_dist = S^{-1}.dx_dist
    El::cholesky::SolveAfter(El::UpperOrLowerNS::UPPER,
                             El::OrientationNS::NORMAL, S, dx_dist);
    // copy dx_dist back to dx
    dx = dx_dist;
  }

  // Compute the search direction (dx, dX, dy, dY) for the predictor and corrector phases.
  //
  // Inputs:
  // - beta, the centering parameter
  // - mu = Tr(X Y) / X.cols
  // - correctorPhase: boolean indicating whether we're in the corrector
  //   phase or predictor phase.
  // Workspace (members of SDPSolver which are modified in-place but not
  // used elsewhere):
  // - Z, R
  // Outputs (members of SDP_Solver which are modified in-place):
  // - dx, dX, dY
  //
  void
  compute_search_direction(const SDP &sdp, const SDP_Solver &solver,
                           const Block_Diagonal_Matrix &minus_XY,
                           const Block_Diagonal_Matrix &X_cholesky,
                           const El::BigFloat &beta, const El::BigFloat &mu,
                           const bool &is_corrector_phase,
                           const El::DistMatrix<El::BigFloat> &S,
                           Primal_Dist_Vector &dx, Block_Diagonal_Matrix &dX,
                           Block_Diagonal_Matrix &dY)
  {
    // R = beta mu I - X Y (predictor phase)
    // R = beta mu I - X Y - dX dY (corrector phase)
    Block_Diagonal_Matrix R(minus_XY);
    if(is_corrector_phase)
      {
        scale_multiply_add(El::BigFloat(-1), dX, dY, El::BigFloat(1), R);
      }
    R.add_diagonal(beta * mu);

    // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
    Block_Diagonal_Matrix Z(solver.X);
    multiply(solver.primal_residues, solver.Y, Z);
    Z -= R;
    cholesky_solve(X_cholesky, Z);
    Z.symmetrize();

    // Solve
    // S_ij dx[j] = -dual_residues[i] - Tr(F_i Z)

    // Write RHS of the equation to dx
    compute_schur_RHS(sdp, solver.dual_residues, Z, dx);
    // Solve for dx in-place
    solve_schur_complement_equation(S, dx);

    // dX = PrimalResidues + \sum_p F_p dx[p]
    constraint_matrix_weighted_sum(sdp, dx, dX);
    dX += solver.primal_residues;

    // dY = Symmetrize(X^{-1} (R - dX Y))
    multiply(dX, solver.Y, dY);
    dY -= R;
    cholesky_solve(X_cholesky, dY);
    dY.symmetrize();
    dY *= El::BigFloat(-1);
  }
}
