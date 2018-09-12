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

void compute_dual_residues(const Block_Info &block_info, const SDP &sdp,
                           const Block_Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Block_Vector &dual_residues);

void compute_primal_residues(const Block_Info &block_info, const SDP &sdp,
                             const Block_Vector &x,
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
SDP_Solver::run(const SDP_Solver_Parameters &parameters,
                const boost::filesystem::path checkpoint_file,
                const Block_Info &block_info, const SDP &sdp,
                const El::Grid &grid, Timers &timers)
{
  SDP_Solver_Terminate_Reason result(
    SDP_Solver_Terminate_Reason::MaxIterationsExceeded);
  auto &solver_timer(timers.add_and_start("Solver runtime"));
  auto &initialize_timer(timers.add_and_start("run.initialize"));

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
  //     0 <= k1,k2 <= block_info.degrees[j] = d_j
  //     0 <= s,r < block_info.dimensions[j] = m_j
  //
  // where j corresponds to b and M^{(s,r)} denotes the (s,r)-th
  // (d_j+1)x(d_j+1) block of M.
  //
  // BilinearPairingsXInv has one block for each block of X.  The
  // dimension of BilinearPairingsXInv.block[b] is (d_j+1)*m_j.  See
  // SDP.h for more information on d_j and m_j.

  Block_Diagonal_Matrix bilinear_pairings_X_inv(
    block_info.bilinear_pairing_block_sizes, block_info.block_indices,
    block_info.schur_block_sizes.size(), grid);

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
          X_block.Height(), bilinear_pairings_X_inv_block->Width(), grid);
        ++bilinear_pairings_X_inv_block;
      }
  }
  print_header();

  std::size_t total_psd_rows(
    std::accumulate(block_info.psd_matrix_block_sizes.begin(),
                    block_info.psd_matrix_block_sizes.end(), size_t(0)));

  initialize_timer.stop();
  boost::timer::cpu_timer checkpoint_timer;
  for(int iteration = 1;; iteration++)
    {
      /// FIXME: This has to use something that is guaranteed to be
      /// the same for all processors.
      if(checkpoint_timer.elapsed().wall
         >= parameters.checkpoint_interval * 1000000000LL)
        {
          save_checkpoint(checkpoint_file);
          checkpoint_timer.start();
        }
      if(solver_timer.elapsed().wall >= parameters.max_runtime * 1000000000LL)
        {
          result = SDP_Solver_Terminate_Reason::MaxRuntimeExceeded;
        }

      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.objectives");
        }
      auto &objectives_timer(timers.add_and_start("run.objectives"));
      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.objectives.primal");
        }
      auto &objectives_primal_timer(
        timers.add_and_start("run.objectives.primal"));
      primal_objective = sdp.objective_const + dot(sdp.primal_objective_c, x);
      objectives_primal_timer.stop();
      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.objectives.dual");
        }
      auto &objectives_dual_timer(timers.add_and_start("run.objectives.dual"));
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

      objectives_dual_timer.stop();
      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.objectives.gap");
        }
      auto &objectives_gap_timer(timers.add_and_start("run.objectives.gap"));
      duality_gap
        = Abs(primal_objective - dual_objective)
          / Max(Abs(primal_objective) + Abs(dual_objective), El::BigFloat(1));
      objectives_gap_timer.stop();

      objectives_timer.stop();

      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(),
                     " run.choleskyDecomposition(X,XCholesky)");
        }
      auto &cholesky_decomposition_X_timer(
        timers.add_and_start("run.choleskyDecomposition(X,XCholesky)"));
      cholesky_decomposition(X, X_cholesky);
      cholesky_decomposition_X_timer.stop();

      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(),
                     " run.choleskyDecomposition(Y,YCholesky)");
        }
      auto &cholesky_decomposition_Y_timer(
        timers.add_and_start("run.choleskyDecomposition(Y,YCholesky)"));
      cholesky_decomposition(Y, Y_cholesky);
      cholesky_decomposition_Y_timer.stop();

      // Compute the bilinear pairings BilinearPairingsXInv and
      // BilinearPairingsY needed for the dualResidues and the Schur
      // complement matrix
      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(),
                     " run.blockTensorInvTransposeCongruenceWithCholesky");
        }
      auto &inv_transpose_congruence_timer(timers.add_and_start(
        "run.blockTensorInvTransposeCongruenceWithCholesky"));
      block_tensor_inv_transpose_congruence_with_cholesky(
        X_cholesky, sdp.bilinear_bases_local, block_info.block_indices,
        bilinear_pairings_workspace, bilinear_pairings_X_inv);
      inv_transpose_congruence_timer.stop();

      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.blockTensorTransposeCongruence");
        }
      auto &transpose_congruence_timer(
        timers.add_and_start("run.blockTensorTransposeCongruence"));
      block_tensor_transpose_congruence(
        Y, sdp.bilinear_bases_local, block_info.block_indices,
        bilinear_pairings_workspace, bilinear_pairings_Y);
      transpose_congruence_timer.stop();

      // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.computeDualResidues");
        }
      auto &dual_residues_timer(
        timers.add_and_start("run.computeDualResidues"));
      compute_dual_residues(block_info, sdp, y, bilinear_pairings_Y,
                            dual_residues);
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
      dual_residues_timer.stop();

      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.computePrimalResidues");
        }
      auto &primal_residues_timer(
        timers.add_and_start("run.computePrimalResidues"));
      // PrimalResidues = \sum_p A_p x[p] - X
      compute_primal_residues(block_info, sdp, x, X, primal_residues);
      primal_error = primal_residues.max_abs();
      primal_residues_timer.stop();

      const bool is_primal_feasible(primal_error
                                    < parameters.primal_error_threshold);
      const bool is_dual_feasible(dual_error
                                  < parameters.dual_error_threshold);
      const bool is_optimal(duality_gap < parameters.duality_gap_threshold);

      if(is_primal_feasible && is_dual_feasible && is_optimal)
        {
          result = SDP_Solver_Terminate_Reason::PrimalDualOptimal;
          break;
        }
      else if(is_primal_feasible && parameters.find_primal_feasible)
        {
          result = SDP_Solver_Terminate_Reason::PrimalFeasible;
          break;
        }
      else if(is_dual_feasible && parameters.find_dual_feasible)
        {
          result = SDP_Solver_Terminate_Reason::DualFeasible;
          break;
        }
      else if(primal_step_length == El::BigFloat(1)
              && parameters.detect_primal_feasible_jump)
        {
          result = SDP_Solver_Terminate_Reason::PrimalFeasibleJumpDetected;
          break;
        }
      else if(dual_step_length == El::BigFloat(1)
              && parameters.detect_dual_feasible_jump)
        {
          result = SDP_Solver_Terminate_Reason::DualFeasibleJumpDetected;
          break;
        }
      else if(iteration > parameters.max_iterations)
        {
          result = SDP_Solver_Terminate_Reason::MaxIterationsExceeded;
          break;
        }

      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.step");
        }
      auto &step_timer(timers.add_and_start("run.step"));
      El::BigFloat mu, beta_predictor, beta_corrector;

      // Search direction: These quantities have the same structure
      // as (x, X, y, Y). They are computed twice each iteration:
      // once in the predictor step, and once in the corrector step.
      Block_Vector dx(x), dy(y);
      Block_Diagonal_Matrix dX(X), dY(Y);
      {
        // SchurComplementCholesky = L', the Cholesky decomposition of the
        // Schur complement matrix S.
        Block_Diagonal_Matrix schur_complement_cholesky(
          block_info.schur_block_sizes, block_info.block_indices,
          block_info.schur_block_sizes.size(), grid);

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
        if(parameters.debug)
          {
            El::Output(El::mpi::Rank(),
                       " run.step.initializeSchurComplementSolver");
          }
        auto &initialize_timer(timers.add_and_start("run.step.initializeSchurComplementSolver"));
        initialize_schur_complement_solver(
          block_info, sdp, bilinear_pairings_X_inv, bilinear_pairings_Y, grid,
          parameters.debug, schur_complement_cholesky, schur_off_diagonal, Q,
          timers);
        initialize_timer.stop();

        // Compute the complementarity mu = Tr(X Y)/X.dim
        if(parameters.debug)
          {
            El::Output(El::mpi::Rank(),
                       " run.step.frobenius_product_symmetric");
          }
        auto &frobenius_timer(timers.add_and_start("run.step.frobenius_product_symmetric"));
        mu = frobenius_product_symmetric(X, Y) / total_psd_rows;
        frobenius_timer.stop();
        if(mu > parameters.max_complementarity)
          {
            result = SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;
            break;
          }

        if(parameters.debug)
          {
            El::Output(El::mpi::Rank(),
                       " run.step.predictor_centering_parameter");
          }
        auto &predictor_timer(timers.add_and_start("run.step.predictor_centering_parameter"));
        // Compute the predictor solution for (dx, dX, dy, dY)
        beta_predictor = predictor_centering_parameter(
          parameters, is_primal_feasible && is_dual_feasible);
        predictor_timer.stop();

        if(parameters.debug)
          {
            El::Output(El::mpi::Rank(),
                       " run.step.computeSearchDirection(betaPredictor)");
          }
        auto &search_predictor_timer(timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));
        compute_search_direction(block_info, sdp, schur_complement_cholesky,
                                 schur_off_diagonal, X_cholesky,
                                 beta_predictor, mu, false, Q, dx, dX, dy, dY);
        search_predictor_timer.stop();

        // Compute the corrector solution for (dx, dX, dy, dY)
        if(parameters.debug)
          {
            El::Output(El::mpi::Rank(),
                       " run.step.corrector_centering_parameter");
          }
        auto &corrector_timer(timers.add_and_start("run.step.corrector_centering_parameter"));
        beta_corrector = corrector_centering_parameter(
          parameters, X, dX, Y, dY, mu, is_primal_feasible && is_dual_feasible,
          total_psd_rows);
        corrector_timer.stop();
        if(parameters.debug)
          {
            El::Output(El::mpi::Rank(),
                       " run.step.computeSearchDirection(betaCorrector)");
          }
        auto &search_corrector_timer(timers.add_and_start("run.step.computeSearchDirection(betaCorrector)"));
        compute_search_direction(block_info, sdp, schur_complement_cholesky,
                                 schur_off_diagonal, X_cholesky,
                                 beta_corrector, mu, true, Q, dx, dX, dy, dY);
        search_corrector_timer.stop();
      }
      // Compute step-lengths that preserve positive definiteness of X, Y
      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.step.stepLength(XCholesky)");
        }
      auto &step_length_X_timer(timers.add_and_start("run.step.stepLength(XCholesky)"));
      primal_step_length
        = step_length(X_cholesky, dX, parameters.step_length_reduction);
      step_length_X_timer.stop();

      if(parameters.debug)
        {
          El::Output(El::mpi::Rank(), " run.step.stepLength(YCholesky)");
        }
      auto &step_length_Y_timer(timers.add_and_start("run.step.stepLength(YCholesky)"));
      dual_step_length
        = step_length(Y_cholesky, dY, parameters.step_length_reduction);
      step_length_Y_timer.stop();

      // If our problem is both dual-feasible and primal-feasible,
      // ensure we're following the true Newton direction.
      if(is_primal_feasible && is_dual_feasible)
        {
          primal_step_length = El::Min(primal_step_length, dual_step_length);
          dual_step_length = primal_step_length;
        }

      print_iteration(iteration, mu, primal_step_length, dual_step_length,
                      beta_corrector, sdp.dual_objective_b.Height(),
                      solver_timer);
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
      step_timer.stop();
    }

  // Never reached
  solver_timer.stop();
  return result;
}
