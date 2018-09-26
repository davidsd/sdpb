#include "../../SDP_Solver.hxx"
#include "../../../Timers.hxx"

// The main solver loop

void print_header();
void print_iteration(
  const int &iteration, const El::BigFloat &mu,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const El::BigFloat &beta_corrector, const size_t &dual_objective_b_height,
  const SDP_Solver &sdp_solver,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time);

void compute_objectives(const SDP &sdp, const Block_Vector &x,
                        const Block_Vector &y, El::BigFloat &primal_objective,
                        El::BigFloat &dual_objective,
                        El::BigFloat &duality_gap, Timers &timers);

void compute_bilinear_pairings(
  const Block_Diagonal_Matrix &X_cholesky, const Block_Diagonal_Matrix &Y,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases,
  std::vector<El::DistMatrix<El::BigFloat>> &workspace,
  Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  Block_Diagonal_Matrix &bilinear_pairings_Y, Timers &timers);

void compute_feasible_and_termination(
  const SDP_Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration, bool &is_primal_and_dual_feasible,
  SDP_Solver_Terminate_Reason &result, bool &terminate_now);

void compute_dual_residues_and_error(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &y,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Vector &dual_residues, El::BigFloat &dual_error, Timers &timers);

void compute_primal_residues_and_error(const Block_Info &block_info,
                                       const SDP &sdp, const Block_Vector &x,
                                       const Block_Diagonal_Matrix &X,
                                       Block_Diagonal_Matrix &primal_residues,
                                       El::BigFloat &primal_error,
                                       Timers &timers);

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
  bilinear_pairings_workspace.reserve(X.blocks.size());
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
  auto last_checkpoint_time(std::chrono::high_resolution_clock::now());
  for(int iteration = 1;; iteration++)
    {
      /// FIXME: This has to use something that is guaranteed to be
      /// the same for all processors.
      if(std::chrono::duration_cast<std::chrono::seconds>(
           std::chrono::high_resolution_clock::now() - last_checkpoint_time)
           .count()
         >= parameters.checkpoint_interval)
        {
          save_checkpoint(checkpoint_file);
          last_checkpoint_time = std::chrono::high_resolution_clock::now();
        }
      if(std::chrono::duration_cast<std::chrono::seconds>(
           std::chrono::high_resolution_clock::now() - solver_timer.start_time)
           .count()
         >= parameters.max_runtime)
        {
          result = SDP_Solver_Terminate_Reason::MaxRuntimeExceeded;
        }

      compute_objectives(sdp, x, y, primal_objective, dual_objective,
                         duality_gap, timers);

      auto &cholesky_decomposition_timer(
        timers.add_and_start("run.choleskyDecomposition"));
      cholesky_decomposition(X, X_cholesky);
      cholesky_decomposition(Y, Y_cholesky);
      cholesky_decomposition_timer.stop();

      compute_bilinear_pairings(
        X_cholesky, Y, sdp.bilinear_bases_local, bilinear_pairings_workspace,
        bilinear_pairings_X_inv, bilinear_pairings_Y, timers);

      compute_dual_residues_and_error(block_info, sdp, y, bilinear_pairings_Y,
                                      dual_residues, dual_error, timers);
      compute_primal_residues_and_error(block_info, sdp, x, X, primal_residues,
                                        primal_error, timers);

      bool is_primal_and_dual_feasible, terminate_now;
      compute_feasible_and_termination(
        parameters, primal_error, dual_error, duality_gap, primal_step_length,
        dual_step_length, iteration, is_primal_and_dual_feasible, result,
        terminate_now);
      if(terminate_now)
        {
          break;
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
        auto &initialize_timer(
          timers.add_and_start("run.step.initializeSchurComplementSolver"));
        initialize_schur_complement_solver(
          block_info, sdp, bilinear_pairings_X_inv, bilinear_pairings_Y, grid,
          schur_complement_cholesky, schur_off_diagonal, Q, timers);
        initialize_timer.stop();

        // Compute the complementarity mu = Tr(X Y)/X.dim
        auto &frobenius_timer(
          timers.add_and_start("run.step.frobenius_product_symmetric"));
        mu = frobenius_product_symmetric(X, Y) / total_psd_rows;
        frobenius_timer.stop();
        if(mu > parameters.max_complementarity)
          {
            result = SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;
            break;
          }

        auto &predictor_timer(
          timers.add_and_start("run.step.predictor_centering_parameter"));
        // Compute the predictor solution for (dx, dX, dy, dY)
        beta_predictor = predictor_centering_parameter(
          parameters, is_primal_and_dual_feasible);
        predictor_timer.stop();

        auto &search_predictor_timer(timers.add_and_start(
          "run.step.computeSearchDirection(betaPredictor)"));
        compute_search_direction(block_info, sdp, schur_complement_cholesky,
                                 schur_off_diagonal, X_cholesky,
                                 beta_predictor, mu, false, Q, dx, dX, dy, dY);
        search_predictor_timer.stop();

        // Compute the corrector solution for (dx, dX, dy, dY)
        auto &corrector_timer(
          timers.add_and_start("run.step.corrector_centering_parameter"));
        beta_corrector = corrector_centering_parameter(
          parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
          total_psd_rows);
        corrector_timer.stop();
        auto &search_corrector_timer(timers.add_and_start(
          "run.step.computeSearchDirection(betaCorrector)"));
        compute_search_direction(block_info, sdp, schur_complement_cholesky,
                                 schur_off_diagonal, X_cholesky,
                                 beta_corrector, mu, true, Q, dx, dX, dy, dY);
        search_corrector_timer.stop();
      }
      // Compute step-lengths that preserve positive definiteness of X, Y
      auto &step_length_X_timer(
        timers.add_and_start("run.step.stepLength(XCholesky)"));
      primal_step_length
        = step_length(X_cholesky, dX, parameters.step_length_reduction);
      step_length_X_timer.stop();

      auto &step_length_Y_timer(
        timers.add_and_start("run.step.stepLength(YCholesky)"));
      dual_step_length
        = step_length(Y_cholesky, dY, parameters.step_length_reduction);
      step_length_Y_timer.stop();

      // If our problem is both dual-feasible and primal-feasible,
      // ensure we're following the true Newton direction.
      if(is_primal_and_dual_feasible)
        {
          primal_step_length = El::Min(primal_step_length, dual_step_length);
          dual_step_length = primal_step_length;
        }

      print_iteration(iteration, mu, primal_step_length, dual_step_length,
                      beta_corrector, sdp.dual_objective_b.Height(), *this,
                      solver_timer.start_time);
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
