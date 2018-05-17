#include "../../SDP_Solver.hxx"
#include "../../../Timers.hxx"

// The main solver loop
// PROFILING : Using the somewhat weird convention to put timings outside the
//             function. Its just simpler to see this way what's timed and
//             what's not

El::BigFloat dot(const Block_Vector &a, const Block_Vector &b);

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &L, const std::vector<Matrix> &Q,
  std::vector<Matrix> &Work, Block_Diagonal_Matrix &Result);

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result);

void block_tensor_transpose_congruence(const Block_Diagonal_Matrix &A,
                                       const std::vector<Matrix> &Q,
                                       std::vector<Matrix> &Work,
                                       Block_Diagonal_Matrix &Result);
void block_tensor_transpose_congruence(
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result);

void compute_dual_residues(const SDP &sdp, const Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Vector &dual_residues);

void compute_dual_residues(const SDP &sdp,
                           const El::DistMatrix<El::BigFloat> &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Block_Vector &dual_residues);

void compute_primal_residues(const SDP &sdp, const Vector &x,
                             const Block_Vector &x_elemental,
                             const Block_Diagonal_Matrix &X,
                             Block_Diagonal_Matrix &primal_residues);

Real predictor_centering_parameter(const SDP_Solver_Parameters &parameters,
                                   const bool is_primal_dual_feasible);
El::BigFloat predictor_centering_parameter_elemental(
  const SDP_Solver_Parameters &parameters, const bool is_primal_dual_feasible);

Real corrector_centering_parameter(const SDP_Solver_Parameters &parameters,
                                   const Block_Diagonal_Matrix &X,
                                   const Block_Diagonal_Matrix &dX,
                                   const Block_Diagonal_Matrix &Y,
                                   const Block_Diagonal_Matrix &dY,
                                   const Real &mu,
                                   const bool is_primal_dual_feasible);

El::BigFloat corrector_centering_parameter(
  const SDP_Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible);

Real step_length(Block_Diagonal_Matrix &MCholesky, Block_Diagonal_Matrix &dM,
                 Block_Diagonal_Matrix &MInvDM,
                 std::vector<Vector> &eigenvalues,
                 std::vector<Vector> &workspace, const Real gamma);
El::BigFloat
step_length(Block_Diagonal_Matrix &MCholesky, Block_Diagonal_Matrix &dM,
            Block_Diagonal_Matrix &MInvDM, const El::BigFloat &gamma);

SDP_Solver_Terminate_Reason
SDP_Solver::run(const boost::filesystem::path checkpoint_file)
{
  Real primal_step_length(0), dual_step_length(0);
  El::BigFloat primal_step_length_elemental(0), dual_step_length_elemental(0);

  // the Cholesky decompositions of X and Y, each lower-triangular
  // BlockDiagonalMatrices with the same block sizes as X and Y
  Block_Diagonal_Matrix X_cholesky(X);
  Block_Diagonal_Matrix Y_cholesky(X);

  // Bilinear pairings needed for computing the Schur complement
  // matrix.  For example,
  //
  //   BilinearPairingsXInv.blocks[b].elt(
  //     (d_j+1) s + k1,
  //     (d_j+1) r + k2
  //   ) = v_{b,k1}^T (X.blocks[b]^{-1})^{(s,r)} v_{b,k2}
  //
  //     0 <= k1,k2 <= sdp.degrees[j] = d_j
  //     0 <= s,r < sdp.dimensions[j] = m_j
  //
  // where j corresponds to b and M^{(s,r)} denotes the (s,r)-th
  // (d_j+1)x(d_j+1) block of M.
  //
  // BilinearPairingsXInv has one block for each block of X.  The
  // dimension of BilinearPairingsXInv.block[b] is (d_j+1)*m_j.  See
  // SDP.h for more information on d_j and m_j.

  Block_Diagonal_Matrix bilinear_pairings_X_Inv(
    sdp.bilinear_pairing_block_dims());

  // BilinearPairingsY is analogous to BilinearPairingsXInv, with
  // X^{-1} -> Y.

  Block_Diagonal_Matrix bilinear_pairings_Y(bilinear_pairings_X_Inv);

  // Additional workspace variables used in step_length()
  Block_Diagonal_Matrix step_matrix_workspace(X);
  std::vector<Matrix> bilinear_pairings_workspace;
  std::vector<El::DistMatrix<El::BigFloat>>
    bilinear_pairings_workspace_elemental;
  std::vector<Vector> eigenvalues_workspace;
  std::vector<Vector> QR_workspace;
  for(unsigned int b = 0; b < sdp.bilinear_bases.size(); b++)
    {
      bilinear_pairings_workspace.emplace_back(
        X.blocks[b].rows, bilinear_pairings_X_Inv.blocks[b].cols);
      bilinear_pairings_workspace_elemental.emplace_back(
        X.blocks_elemental[b].Height(),
        bilinear_pairings_X_Inv.blocks_elemental[b].Width());
      eigenvalues_workspace.emplace_back(X.blocks[b].rows);
      QR_workspace.emplace_back(3 * X.blocks[b].rows - 1);
    }

  print_header();

  for(int iteration = 1;; iteration++)
    {
      if(timers["Last checkpoint"].elapsed().wall
         >= parameters.checkpoint_interval * 1000000000LL)
        {
          save_checkpoint(checkpoint_file);
        }
      if(timers["Solver runtime"].elapsed().wall
         >= parameters.max_runtime * 1000000000LL)
        {
          return SDP_Solver_Terminate_Reason::MaxRuntimeExceeded;
        }

      timers["run.objectives"].resume();
      primal_objective
        = sdp.objective_const + dot_product(sdp.primal_objective_c, x);
      dual_objective
        = sdp.objective_const + dot_product(sdp.dual_objective_b, y);
      duality_gap
        = abs(primal_objective - dual_objective)
          / max(Real(abs(primal_objective) + abs(dual_objective)), Real(1));

      primal_objective_elemental
        = sdp.objective_const_elemental
          + dot(sdp.primal_objective_c_elemental, x_elemental);
      dual_objective_elemental
        = sdp.objective_const_elemental
          + El::Dotu(sdp.dual_objective_b_elemental, y_elemental);
      duality_gap_elemental
        = Abs(primal_objective_elemental - dual_objective_elemental)
          / Max(Abs(primal_objective_elemental)
                  + Abs(dual_objective_elemental),
                El::BigFloat(1));

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

      block_tensor_inv_transpose_congruence_with_cholesky(
        X_cholesky, sdp.bilinear_bases_elemental_local,
        bilinear_pairings_workspace_elemental, bilinear_pairings_X_Inv);
      timers["run.blockTensorInvTransposeCongruenceWithCholesky"].stop();

      timers["run.blockTensorTransposeCongruence"].resume();
      block_tensor_transpose_congruence(Y, sdp.bilinear_bases,
                                        bilinear_pairings_workspace,
                                        bilinear_pairings_Y);

      block_tensor_transpose_congruence(Y, sdp.bilinear_bases_elemental_dist,
                                        bilinear_pairings_workspace_elemental,
                                        bilinear_pairings_Y);
      timers["run.blockTensorTransposeCongruence"].stop();

      // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
      timers["run.computeDualResidues"].resume();
      compute_dual_residues(sdp, y, bilinear_pairings_Y, dual_residues);
      dual_error = max_abs_vector(dual_residues);

      compute_dual_residues(sdp, y_elemental, bilinear_pairings_Y,
                            dual_residues_elemental);
      dual_error_elemental = 0;
      for(auto &block : dual_residues_elemental.blocks)
        {
          dual_error_elemental = Max(dual_error_elemental, El::MaxAbs(block));
        }
      timers["run.computeDualResidues"].stop();

      timers["run.computePrimalResidues"].resume();
      // PrimalResidues = \sum_p A_p x[p] - X
      compute_primal_residues(sdp, x, x_elemental, X, primal_residues);
      primal_error_elemental = primal_residues.max_abs();
      timers["run.computePrimalResidues"].stop();

      const bool is_primal_feasible(
        primal_error_elemental < parameters.primal_error_threshold_elemental);
      const bool is_dual_feasible(dual_error_elemental
                                  < parameters.dual_error_threshold_elemental);
      const bool is_optimal(duality_gap_elemental
                            < parameters.duality_gap_threshold_elemental);

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

      Real mu, beta_predictor, beta_corrector;
      El::BigFloat mu_elemental, beta_predictor_elemental,
        beta_corrector_elemental;

      // Search direction: These quantities have the same structure
      // as (x, X, y, Y). They are computed twice each iteration:
      // once in the predictor step, and once in the corrector step.
      Vector dx(x), dy(y);
      Block_Vector dx_elemental(x_elemental);
      El::DistMatrix<El::BigFloat> dy_elemental(y_elemental);
      Block_Diagonal_Matrix dX(X), dY(Y);
      {
        // FIXME: It may be expensive to create these objects for each
        // iteration.

        // SchurComplementCholesky = L', the Cholesky decomposition of the
        // Schur complement matrix S.
        Block_Diagonal_Matrix schur_complement_cholesky(
          sdp.schur_block_dims());

        // SchurOffDiagonal = L'^{-1} FreeVarMatrix, needed in solving the
        // Schur complement equation.
        Matrix schur_off_diagonal;
        Block_Matrix schur_off_diagonal_elemental;

        // Q = B' L'^{-T} L'^{-1} B' - {{0, 0}, {0, 1}}, where B' =
        // (FreeVarMatrix U).  Q is needed in the factorization of the Schur
        // complement equation.  Q has dimension N'xN', where
        //
        //   N' = cols(B) + cols(U) = N + cols(U)
        //
        // where N is the dimension of the dual objective function.  Note
        // that N' could change with each iteration.
        Matrix Q(sdp.free_var_matrix.cols, sdp.free_var_matrix.cols);
        El::DistMatrix<El::BigFloat> Q_elemental(sdp.free_var_matrix.cols,
                                                 sdp.free_var_matrix.cols);

        // Compute SchurComplement and prepare to solve the Schur
        // complement equation for dx, dy
        timers["run.initializeSchurComplementSolver"].resume();
        initialize_schur_complement_solver(
          bilinear_pairings_X_Inv, bilinear_pairings_Y, sdp.schur_block_dims(),
          schur_complement_cholesky, schur_off_diagonal,
          schur_off_diagonal_elemental, Q, Q_elemental);
        timers["run.initializeSchurComplementSolver"].stop();

        // Compute the complementarity mu = Tr(X Y)/X.dim
        mu = frobenius_product_symmetric(X, Y) / X.dim;
        if(mu > parameters.max_complementarity)
          return SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;

        mu_elemental = frobenius_product_symmetric_elemental(X, Y) / X.dim;
        if(mu_elemental > parameters.max_complementarity_elemental)
          {
            return SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;
          }

        // Compute the predictor solution for (dx, dX, dy, dY)
        beta_predictor = predictor_centering_parameter(
          parameters, is_primal_feasible && is_dual_feasible);

        beta_predictor_elemental = predictor_centering_parameter_elemental(
          parameters, is_primal_feasible && is_dual_feasible);

        timers["run.computeSearchDirection(betaPredictor)"].resume();
        compute_search_direction(schur_complement_cholesky, schur_off_diagonal,
                                 schur_off_diagonal_elemental, X_cholesky,
                                 beta_predictor, beta_predictor_elemental, mu,
                                 mu_elemental, false, Q, Q_elemental, dx,
                                 dx_elemental, dX, dy, dy_elemental, dY);
        timers["run.computeSearchDirection(betaPredictor)"].stop();

        // Compute the corrector solution for (dx, dX, dy, dY)
        beta_corrector_elemental = corrector_centering_parameter(
          parameters, X, dX, Y, dY, mu_elemental,
          is_primal_feasible && is_dual_feasible);
        timers["run.computeSearchDirection(betaCorrector)"].resume();
        compute_search_direction(schur_complement_cholesky, schur_off_diagonal,
                                 schur_off_diagonal_elemental, X_cholesky,
                                 beta_corrector, beta_corrector_elemental, mu,
                                 mu_elemental, true, Q, Q_elemental, dx,
                                 dx_elemental, dX, dy, dy_elemental, dY);
        timers["run.computeSearchDirection(betaCorrector)"].stop();
      }
      // Compute step-lengths that preserve positive definiteness of X, Y
      timers["run.stepLength(XCholesky)"].resume();
      primal_step_length_elemental
        = step_length(X_cholesky, dX, step_matrix_workspace,
                      parameters.step_length_reduction_elemental);

      timers["run.stepLength(XCholesky)"].stop();
      timers["run.stepLength(YCholesky)"].resume();
      dual_step_length_elemental
        = step_length(Y_cholesky, dY, step_matrix_workspace,
                      parameters.step_length_reduction_elemental);
      timers["run.stepLength(YCholesky)"].stop();

      // If our problem is both dual-feasible and primal-feasible,
      // ensure we're following the true Newton direction.
      if(is_primal_feasible && is_dual_feasible)
        {
          primal_step_length_elemental
            = min(primal_step_length_elemental, dual_step_length_elemental);
          dual_step_length_elemental = primal_step_length_elemental;
        }

      print_iteration(iteration, mu_elemental, primal_step_length_elemental,
                      dual_step_length_elemental, beta_corrector_elemental);
      // Update the primal point (x, X) += primalStepLength*(dx, dX)
      for(size_t block = 0; block < x_elemental.blocks.size(); ++block)
        {
          El::Axpy(primal_step_length_elemental, dx_elemental.blocks[block],
                   x_elemental.blocks[block]);
        }
      dX *= primal_step_length_elemental;

      X += dX;

      // Update the dual point (y, Y) += dualStepLength*(dy, dY)
      El::Axpy(dual_step_length_elemental, dy_elemental, y_elemental);
      dY *= dual_step_length_elemental;

      Y += dY;
    }

  // Never reached
  return SDP_Solver_Terminate_Reason::MaxIterationsExceeded;
}
