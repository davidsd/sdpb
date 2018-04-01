#include "../../SDP_Solver.hxx"
#include "../../../Timers.hxx"

// The main solver loop
// PROFILING : Using the somewhat weird convention to put timings outside the
//             function. Its just simpler to see this way what's timed and
//             what's not

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &L, const std::vector<Matrix> &Q,
  std::vector<Matrix> &Work, Block_Diagonal_Matrix &Result);

void block_tensor_transpose_congruence(const Block_Diagonal_Matrix &A,
                                       const std::vector<Matrix> &Q,
                                       std::vector<Matrix> &Work,
                                       Block_Diagonal_Matrix &Result);

void compute_dual_residues(const SDP &sdp, const Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Vector &dual_residues);

void compute_primal_residues(const SDP &sdp, const Vector x,
                             const Block_Diagonal_Matrix &X,
                             Block_Diagonal_Matrix &primal_residues);

Real predictor_centering_parameter(const SDP_Solver_Parameters &parameters,
                                   const bool is_primal_dual_feasible);

Real corrector_centering_parameter(const SDP_Solver_Parameters &parameters,
                                   const Block_Diagonal_Matrix &X,
                                   const Block_Diagonal_Matrix &dX,
                                   const Block_Diagonal_Matrix &Y,
                                   const Block_Diagonal_Matrix &dY,
                                   const Real &mu,
                                   const bool is_primal_dual_feasible);

Real step_length(Block_Diagonal_Matrix &MCholesky, Block_Diagonal_Matrix &dM,
                 Block_Diagonal_Matrix &MInvDM,
                 std::vector<Vector> &eigenvalues,
                 std::vector<Vector> &workspace, const Real gamma);

SDP_Solver_Terminate_Reason
SDP_Solver::run(const boost::filesystem::path checkpoint_file)
{
  Real primal_step_length;
  Real dual_step_length;

  print_header();

  for(int iteration = 1;; iteration++)
    {
      if(timers["Last checkpoint"].elapsed().wall
         >= parameters.checkpoint_interval * 1000000000LL)
        save_checkpoint(checkpoint_file);
      if(timers["Solver runtime"].elapsed().wall
         >= parameters.max_runtime * 1000000000LL)
        return SDP_Solver_Terminate_Reason::MaxRuntimeExceeded;

      timers["run.objectives"].resume();
      primal_objective
        = sdp.objective_const + dot_product(sdp.primal_objective, x);
      dual_objective
        = sdp.objective_const + dot_product(sdp.dual_objective, y);
      duality_gap
        = abs(primal_objective - dual_objective)
          / max(Real(abs(primal_objective) + abs(dual_objective)), Real(1));
      timers["run.objectives"].stop();

      timers["run.choleskyDecomposition(X,XCholesky)"].resume();
      cholesky_decomposition(X, X_cholesky);
      timers["run.choleskyDecomposition(X,XCholesky)"].stop();

      timers["run.choleskyDecomposition(Y,YCholesky)"].resume();
      cholesky_decomposition(Y, Y_cholesky);
      timers["run.choleskyDecomposition(Y,YCholesky)"].stop();

      // Compute the bilinear pairings BilinearPairingsXInv and
      // BilinearPairingsY needed for the dualResidues and the Schur
      // complement matrix
      timers["run.blockTensorInvTransposeCongruenceWithCholesky"].resume();
      block_tensor_inv_transpose_congruence_with_cholesky(
        X_cholesky, sdp.bilinear_bases, bilinear_pairings_workspace,
        bilinear_pairings_X_Inv);
      timers["run.blockTensorInvTransposeCongruenceWithCholesky"].stop();

      timers["run.blockTensorTransposeCongruence"].resume();
      block_tensor_transpose_congruence(Y, sdp.bilinear_bases,
                                        bilinear_pairings_workspace,
                                        bilinear_pairings_Y);
      timers["run.blockTensorTransposeCongruence"].stop();

      // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
      timers["run.computeDualResidues"].resume();
      compute_dual_residues(sdp, y, bilinear_pairings_Y, dual_residues);
      dual_error = max_abs_vector(dual_residues);
      timers["run.computeDualResidues"].stop();

      timers["run.computePrimalResidues"].resume();
      // PrimalResidues = \sum_p A_p x[p] - X
      compute_primal_residues(sdp, x, X, primal_residues);
      primal_error = primal_residues.max_abs();
      timers["run.computePrimalResidues"].stop();

      const bool is_primal_feasible
        = primal_error < parameters.primal_error_threshold;
      const bool is_dual_feasible
        = dual_error < parameters.dual_error_threshold;
      const bool is_optimal = duality_gap < parameters.duality_gap_threshold;

      if(is_primal_feasible && is_dual_feasible && is_optimal)
        return SDP_Solver_Terminate_Reason::PrimalDualOptimal;
      else if(is_primal_feasible && parameters.find_primal_feasible)
        return SDP_Solver_Terminate_Reason::PrimalFeasible;
      else if(is_dual_feasible && parameters.find_dual_feasible)
        return SDP_Solver_Terminate_Reason::DualFeasible;
      else if(primal_step_length == 1
              && parameters.detect_primal_feasible_jump)
        return SDP_Solver_Terminate_Reason::PrimalFeasibleJumpDetected;
      else if(dual_step_length == 1 && parameters.detect_dual_feasible_jump)
        return SDP_Solver_Terminate_Reason::DualFeasibleJumpDetected;
      else if(iteration > parameters.max_iterations)
        return SDP_Solver_Terminate_Reason::MaxIterationsExceeded;

      // Compute SchurComplement and prepare to solve the Schur
      // complement equation for dx, dy
      timers["run.initializeSchurComplementSolver"].resume();
      initialize_schur_complement_solver(bilinear_pairings_X_Inv,
                                         bilinear_pairings_Y);
      timers["run.initializeSchurComplementSolver"].stop();

      // Compute the complementarity mu = Tr(X Y)/X.dim
      Real mu = frobenius_product_symmetric(X, Y) / X.dim;
      if(mu > parameters.max_complementarity)
        return SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;

      // Compute the predictor solution for (dx, dX, dy, dY)
      Real beta_predictor = predictor_centering_parameter(
        parameters, is_primal_feasible && is_dual_feasible);
      timers["run.computeSearchDirection(betaPredictor)"].resume();
      compute_search_direction(beta_predictor, mu, false);
      timers["run.computeSearchDirection(betaPredictor)"].stop();

      // Compute the corrector solution for (dx, dX, dy, dY)
      Real beta_corrector = corrector_centering_parameter(
        parameters, X, dX, Y, dY, mu, is_primal_feasible && is_dual_feasible);
      timers["run.computeSearchDirection(betaCorrector)"].resume();
      compute_search_direction(beta_corrector, mu, true);
      timers["run.computeSearchDirection(betaCorrector)"].stop();

      // Compute step-lengths that preserve positive definiteness of X, Y
      timers["run.stepLength(XCholesky)"].resume();
      primal_step_length = step_length(X_cholesky, dX, step_matrix_workspace,
                                       eigenvalues_workspace, QR_workspace,
                                       parameters.step_length_reduction);
      timers["run.stepLength(XCholesky)"].stop();
      timers["run.stepLength(YCholesky)"].resume();
      dual_step_length = step_length(Y_cholesky, dY, step_matrix_workspace,
                                     eigenvalues_workspace, QR_workspace,
                                     parameters.step_length_reduction);
      timers["run.stepLength(YCholesky)"].stop();

      // If our problem is both dual-feasible and primal-feasible,
      // ensure we're following the true Newton direction.
      if(is_primal_feasible && is_dual_feasible)
        {
          primal_step_length = min(primal_step_length, dual_step_length);
          dual_step_length = primal_step_length;
        }

      print_iteration(iteration, mu, primal_step_length, dual_step_length,
                      beta_corrector);

      // Update the primal point (x, X) += primalStepLength*(dx, dX)
      add_scaled_vector(x, primal_step_length, dx);
      dX *= primal_step_length;
      X += dX;

      // Update the dual point (y, Y) += dualStepLength*(dy, dY)
      add_scaled_vector(y, dual_step_length, dy);
      dY *= dual_step_length;
      Y += dY;
    }

  // Never reached
  return SDP_Solver_Terminate_Reason::MaxIterationsExceeded;
}
