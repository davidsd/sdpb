#include "../../SDP_Solver.hxx"

// The main solver loop

void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);

void print_header(const Verbosity &verbosity);
void print_iteration(
  const int &iteration, const El::BigFloat &mu,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const El::BigFloat &beta_corrector, const SDP_Solver &sdp_solver,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  const Verbosity &verbosity);

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
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  SDP_Solver_Terminate_Reason &terminate_reason, bool &terminate_now);

void compute_dual_residues_and_error(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &y,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Vector &dual_residues, El::BigFloat &dual_error, Timers &timers);

void compute_primal_residues_and_error_P_Ax_X(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &x,
  const Block_Diagonal_Matrix &X, Block_Diagonal_Matrix &primal_residues,
  El::BigFloat &primal_error_P, Timers &timers);

void compute_primal_residues_and_error_p_b_Bx(const Block_Info &block_info,
                                              const SDP &sdp,
                                              const Block_Vector &x,
                                              Block_Vector &primal_residue_p,
                                              El::BigFloat &primal_error_p);

SDP_Solver_Terminate_Reason
SDP_Solver::run(const SDP_Solver_Parameters &parameters,
                const Block_Info &block_info, const SDP &sdp,
                const El::Grid &grid, Timers &timers)
{
  SDP_Solver_Terminate_Reason terminate_reason(
    SDP_Solver_Terminate_Reason::MaxIterationsExceeded);
  auto &solver_timer(timers.add_and_start("Solver runtime"));
  auto &initialize_timer(timers.add_and_start("run.initialize"));

  El::BigFloat primal_step_length(0), dual_step_length(0);

  Block_Diagonal_Matrix X_cholesky(X), Y_cholesky(X);

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
    block_info.bilinear_pairing_block_sizes(), block_info.block_indices,
    block_info.num_points.size(), grid);

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
  print_header(parameters.verbosity);

  auto psd_sizes(block_info.psd_matrix_block_sizes());
  std::size_t total_psd_rows(
    std::accumulate(psd_sizes.begin(), psd_sizes.end(), size_t(0)));

  initialize_timer.stop();
  auto last_checkpoint_time(std::chrono::high_resolution_clock::now());
  for(size_t iteration = 1;; ++iteration)
    {
      El::byte checkpoint_now(
        std::chrono::duration_cast<std::chrono::seconds>(
          std::chrono::high_resolution_clock::now() - last_checkpoint_time)
          .count()
        >= parameters.checkpoint_interval);
      // Time varies between cores, so follow the decision of the root.
      El::mpi::Broadcast(checkpoint_now, 0, El::mpi::COMM_WORLD);
      if(checkpoint_now == true)
        {
          save_checkpoint(parameters);
          last_checkpoint_time = std::chrono::high_resolution_clock::now();
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
      compute_primal_residues_and_error_P_Ax_X(
        block_info, sdp, x, X, primal_residues, primal_error_P, timers);

      // use y to set the sizes of primal_residue_p.  The data is
      // overwritten in compute_primal_residues_and_error_p.
      Block_Vector primal_residue_p(y);
      compute_primal_residues_and_error_p_b_Bx(
        block_info, sdp, x, primal_residue_p, primal_error_p);

      bool terminate_now, is_primal_and_dual_feasible;
      compute_feasible_and_termination(
        parameters, primal_error(), dual_error, duality_gap,
        primal_step_length, dual_step_length, iteration,
        solver_timer.start_time, is_primal_and_dual_feasible, terminate_reason,
        terminate_now);
      if(terminate_now)
        {
          break;
        }

      El::BigFloat mu, beta_corrector;
      step(parameters, total_psd_rows, is_primal_and_dual_feasible, block_info,
           sdp, grid, X_cholesky, Y_cholesky, bilinear_pairings_X_inv,
           bilinear_pairings_Y, primal_residue_p, mu, beta_corrector,
           primal_step_length, dual_step_length, terminate_now, timers);
      if(terminate_now)
        {
          terminate_reason
            = SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;
          break;
        }
      print_iteration(iteration, mu, primal_step_length, dual_step_length,
                      beta_corrector, *this, solver_timer.start_time,
                      parameters.verbosity);
    }
  solver_timer.stop();
  return terminate_reason;
}
