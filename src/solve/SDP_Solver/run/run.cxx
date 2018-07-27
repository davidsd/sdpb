#include "../../SDP_Solver.hxx"
#include "../../../Timers.hxx"

// The main solver loop
// PROFILING : Using the somewhat weird convention to put timings outside the
//             function. Its just simpler to see this way what's timed and
//             what's not

El::BigFloat dot(const Block_Vector &a, const Block_Vector &b);

void block_tensor_inv_transpose_congruence_with_cholesky(
  const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  const std::vector<size_t> &block_indices,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result);

void block_tensor_transpose_congruence(
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  const std::vector<size_t> &block_indices,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &result);

void compute_dual_residues(const SDP &sdp, const Block_Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Block_Vector &dual_residues);

void compute_primal_residues(const SDP &sdp, const Block_Vector &x,
                             const Block_Diagonal_Matrix &X,
                             Block_Diagonal_Matrix &primal_residues);

El::BigFloat
predictor_centering_parameter(const SDP_Solver_Parameters &parameters,
                              const bool is_primal_dual_feasible);

El::BigFloat corrector_centering_parameter(
  const SDP_Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_num_rows);

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma);

SDP_Solver_Terminate_Reason
SDP_Solver::run(const boost::filesystem::path checkpoint_file,
                const bool &debug)
{
  timers["run.initialize"].resume();
  El::BigFloat primal_step_length(0), dual_step_length(0);

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

  Block_Diagonal_Matrix bilinear_pairings_X_inv(
    sdp.bilinear_pairing_block_sizes, sdp.block_indices,
    sdp.schur_block_sizes.size(), sdp.grid);

  // BilinearPairingsY is analogous to BilinearPairingsXInv, with
  // X^{-1} -> Y.

  Block_Diagonal_Matrix bilinear_pairings_Y(bilinear_pairings_X_inv);

  // Additional workspace variables used in step_length()
  std::vector<El::DistMatrix<El::BigFloat>> bilinear_pairings_workspace;
  {
    auto bilinear_pairings_X_inv_block(bilinear_pairings_X_inv.blocks.begin());
    for(auto &X_block : X.blocks)
      {
        bilinear_pairings_workspace.emplace_back(
          X_block.Height(), bilinear_pairings_X_inv_block->Width(), sdp.grid);
        ++bilinear_pairings_X_inv_block;
      }
  }
  print_header();

  timers["run.initialize"].stop();
  for(int iteration = 1;; iteration++)
    {
      /// FIXME: This has to use something that is guaranteed to be
      /// the same for all processors.
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

      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.objectives");
        }
      timers["run.objectives"].resume();
      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.objectives.primal");
        }
      timers["run.objectives.primal"].resume();
      primal_objective = sdp.objective_const + dot(sdp.primal_objective_c, x);
      timers["run.objectives.primal"].stop();
      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.objectives.dual");
        }
      timers["run.objectives.dual"].resume();
      // dual_objective_b is duplicated amongst the processors.  y is
      // duplicated amongst the blocks, but it is possible for some
      // processors to have no blocks.  In principle, we only need to
      // compute the dot product on the first block, but then we would
      // have to make sure that we compute that product over all
      // processors that own that block.
      if(!y.blocks.empty())
        {
          dual_objective = sdp.objective_const
                           + El::Dotu(sdp.dual_objective_b, y.blocks.front());
        }
      El::mpi::Broadcast(&dual_objective, 1, 0, El::mpi::COMM_WORLD);

      timers["run.objectives.dual"].stop();
      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.objectives.gap");
        }
      timers["run.objectives.gap"].resume();
      duality_gap
        = Abs(primal_objective - dual_objective)
          / Max(Abs(primal_objective) + Abs(dual_objective), El::BigFloat(1));
      timers["run.objectives.gap"].stop();

      timers["run.objectives"].stop();

      if(debug)
        {
          El::Output(El::mpi::Rank(),
                     "run.choleskyDecomposition(X,XCholesky)");
        }
      timers["run.choleskyDecomposition(X,XCholesky)"].resume();
      cholesky_decomposition(X, X_cholesky);
      timers["run.choleskyDecomposition(X,XCholesky)"].stop();

      if(debug)
        {
          El::Output(El::mpi::Rank(),
                     "run.choleskyDecomposition(Y,YCholesky)");
        }
      timers["run.choleskyDecomposition(Y,YCholesky)"].resume();
      cholesky_decomposition(Y, Y_cholesky);
      timers["run.choleskyDecomposition(Y,YCholesky)"].stop();

      // Compute the bilinear pairings BilinearPairingsXInv and
      // BilinearPairingsY needed for the dualResidues and the Schur
      // complement matrix
      if(debug)
        {
          El::Output(El::mpi::Rank(),
                     "run.blockTensorInvTransposeCongruenceWithCholesky");
        }
      timers["run.blockTensorInvTransposeCongruenceWithCholesky"].resume();
      block_tensor_inv_transpose_congruence_with_cholesky(
        X_cholesky, sdp.bilinear_bases_local, sdp.block_indices,
        bilinear_pairings_workspace, bilinear_pairings_X_inv);

      timers["run.blockTensorInvTransposeCongruenceWithCholesky"].stop();

      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.blockTensorTransposeCongruence");
        }
      timers["run.blockTensorTransposeCongruence"].resume();
      block_tensor_transpose_congruence(
        Y, sdp.bilinear_bases_local, sdp.block_indices,
        bilinear_pairings_workspace, bilinear_pairings_Y);
      timers["run.blockTensorTransposeCongruence"].stop();

      // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.computeDualResidues");
        }
      timers["run.computeDualResidues"].resume();
      compute_dual_residues(sdp, y, bilinear_pairings_Y, dual_residues);
      {
        dual_error = 0;
        El::BigFloat local_max(0);
        for(auto &block : dual_residues.blocks)
          {
            local_max = Max(local_max, El::MaxAbs(block));
          }
        dual_error
          = El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);
      }
      timers["run.computeDualResidues"].stop();

      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.computePrimalResidues");
        }
      timers["run.computePrimalResidues"].resume();
      // PrimalResidues = \sum_p A_p x[p] - X
      compute_primal_residues(sdp, x, X, primal_residues);
      primal_error = primal_residues.max_abs();
      timers["run.computePrimalResidues"].stop();

      const bool is_primal_feasible(primal_error
                                    < parameters.primal_error_threshold);
      const bool is_dual_feasible(dual_error
                                  < parameters.dual_error_threshold);
      const bool is_optimal(duality_gap < parameters.duality_gap_threshold);

      if(is_primal_feasible && is_dual_feasible && is_optimal)
        return SDP_Solver_Terminate_Reason::PrimalDualOptimal;
      else if(is_primal_feasible && parameters.find_primal_feasible)
        return SDP_Solver_Terminate_Reason::PrimalFeasible;
      else if(is_dual_feasible && parameters.find_dual_feasible)
        return SDP_Solver_Terminate_Reason::DualFeasible;
      else if(primal_step_length == El::BigFloat(1)
              && parameters.detect_primal_feasible_jump)
        return SDP_Solver_Terminate_Reason::PrimalFeasibleJumpDetected;
      else if(dual_step_length == El::BigFloat(1)
              && parameters.detect_dual_feasible_jump)
        return SDP_Solver_Terminate_Reason::DualFeasibleJumpDetected;
      else if(iteration > parameters.max_iterations)
        return SDP_Solver_Terminate_Reason::MaxIterationsExceeded;

      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.step");
        }
      timers["run.step"].resume();
      El::BigFloat mu, beta_predictor, beta_corrector;

      // Search direction: These quantities have the same structure
      // as (x, X, y, Y). They are computed twice each iteration:
      // once in the predictor step, and once in the corrector step.
      Block_Vector dx(x), dy(y);
      Block_Diagonal_Matrix dX(X), dY(Y);
      {
        // FIXME: It may be expensive to create these objects for each
        // iteration.

        // SchurComplementCholesky = L', the Cholesky decomposition of the
        // Schur complement matrix S.
        Block_Diagonal_Matrix schur_complement_cholesky(
          sdp.schur_block_sizes, sdp.block_indices,
          sdp.schur_block_sizes.size(), sdp.grid);

        // SchurOffDiagonal = L'^{-1} FreeVarMatrix, needed in solving the
        // Schur complement equation.
        Block_Matrix schur_off_diagonal;

        // Q = B' L'^{-T} L'^{-1} B' - {{0, 0}, {0, 1}}, where B' =
        // (FreeVarMatrix U).  Q is needed in the factorization of the Schur
        // complement equation.  Q has dimension N'xN', where
        //
        //   N' = cols(B) + cols(U) = N + cols(U)
        //
        // where N is the dimension of the dual objective function.  Note
        // that N' could change with each iteration.
        El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
                                       sdp.dual_objective_b.Height());

        // Compute SchurComplement and prepare to solve the Schur
        // complement equation for dx, dy
        if(debug)
          {
            El::Output(El::mpi::Rank(),
                       "run.step.initializeSchurComplementSolver");
          }
        timers["run.step.initializeSchurComplementSolver"].resume();
        initialize_schur_complement_solver(
          bilinear_pairings_X_inv, bilinear_pairings_Y, sdp.schur_block_sizes,
          sdp.block_indices, sdp.schur_block_sizes.size(), sdp.grid,
          sdp.schur_offsets(), debug, schur_complement_cholesky,
          schur_off_diagonal, Q);
        timers["run.step.initializeSchurComplementSolver"].stop();

        // Compute the complementarity mu = Tr(X Y)/X.dim
        if(debug)
          {
            El::Output(El::mpi::Rank(),
                       "run.step.frobenius_product_symmetric");
          }
        timers["run.step.frobenius_product_symmetric"].resume();
        mu = frobenius_product_symmetric(X, Y) / sdp.total_psd_rows();
        timers["run.step.frobenius_product_symmetric"].stop();
        if(mu > parameters.max_complementarity)
          {
            return SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;
          }

        if(debug)
          {
            El::Output(El::mpi::Rank(),
                       "run.step.predictor_centering_parameter");
          }
        timers["run.step.predictor_centering_parameter"].resume();
        // Compute the predictor solution for (dx, dX, dy, dY)
        beta_predictor = predictor_centering_parameter(
          parameters, is_primal_feasible && is_dual_feasible);
        timers["run.step.predictor_centering_parameter"].stop();

        if(debug)
          {
            El::Output(El::mpi::Rank(),
                       "run.step.computeSearchDirection(betaPredictor)");
          }
        timers["run.step.computeSearchDirection(betaPredictor)"].resume();
        compute_search_direction(schur_complement_cholesky, schur_off_diagonal,
                                 X_cholesky, beta_predictor, mu, false, Q, dx,
                                 dX, dy, dY);
        timers["run.step.computeSearchDirection(betaPredictor)"].stop();

        // Compute the corrector solution for (dx, dX, dy, dY)
        if(debug)
          {
            El::Output(El::mpi::Rank(),
                       "run.step.corrector_centering_parameter");
          }
        timers["run.step.corrector_centering_parameter"].resume();
        beta_corrector = corrector_centering_parameter(
          parameters, X, dX, Y, dY, mu, is_primal_feasible && is_dual_feasible,
          sdp.total_psd_rows());
        timers["run.step.corrector_centering_parameter"].stop();
        if(debug)
          {
            El::Output(El::mpi::Rank(),
                       "run.step.computeSearchDirection(betaCorrector)");
          }
        timers["run.step.computeSearchDirection(betaCorrector)"].resume();
        compute_search_direction(schur_complement_cholesky, schur_off_diagonal,
                                 X_cholesky, beta_corrector, mu, true, Q, dx,
                                 dX, dy, dY);
        timers["run.step.computeSearchDirection(betaCorrector)"].stop();
      }
      // Compute step-lengths that preserve positive definiteness of X, Y
      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.step.stepLength(XCholesky)");
        }
      timers["run.step.stepLength(XCholesky)"].resume();
      primal_step_length
        = step_length(X_cholesky, dX, parameters.step_length_reduction);

      timers["run.step.stepLength(XCholesky)"].stop();

      if(debug)
        {
          El::Output(El::mpi::Rank(), "run.step.stepLength(YCholesky)");
        }
      timers["run.step.stepLength(YCholesky)"].resume();
      dual_step_length
        = step_length(Y_cholesky, dY, parameters.step_length_reduction);
      timers["run.step.stepLength(YCholesky)"].stop();

      // If our problem is both dual-feasible and primal-feasible,
      // ensure we're following the true Newton direction.
      if(is_primal_feasible && is_dual_feasible)
        {
          primal_step_length = El::Min(primal_step_length, dual_step_length);
          dual_step_length = primal_step_length;
        }

      print_iteration(iteration, mu, primal_step_length, dual_step_length,
                      beta_corrector);
      // Update the primal point (x, X) += primalStepLength*(dx, dX)
      for(size_t block = 0; block < x.blocks.size(); ++block)
        {
          El::Axpy(primal_step_length, dx.blocks[block], x.blocks[block]);
        }
      dX *= primal_step_length;

      X += dX;

      // Update the dual point (y, Y) += dualStepLength*(dy, dY)
      for(size_t block = 0; block < dy.blocks.size(); ++block)
        {
          El::Axpy(dual_step_length, dy.blocks[block], y.blocks[block]);
        }
      dY *= dual_step_length;

      Y += dY;
      timers["run.step"].stop();
    }

  // Never reached
  return SDP_Solver_Terminate_Reason::MaxIterationsExceeded;
}
