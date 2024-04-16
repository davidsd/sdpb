#include "sdp_solve/SDP_Solver.hxx"
#include "bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "bigint_syrk/initialize_bigint_syrk_context.hxx"
#include <boost/date_time/posix_time/posix_time.hpp>

// The main solver loop

namespace fs = std::filesystem;

void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L,
                            const Block_Info &block_info,
                            const std::string &name);

void print_header(const Verbosity &verbosity);
void print_iteration(
  const fs::path &iterations_json_path, const int &iteration,
  const El::BigFloat &mu, const El::BigFloat &primal_step_length,
  const El::BigFloat &dual_step_length, const El::BigFloat &beta_corrector,
  const SDP_Solver &sdp_solver,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  const Verbosity &verbosity);

void compute_objectives(const SDP &sdp, const Block_Vector &x,
                        const Block_Vector &y, El::BigFloat &primal_objective,
                        El::BigFloat &dual_objective,
                        El::BigFloat &duality_gap, Timers &timers);

void compute_bilinear_pairings(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_X_inv,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_Y,
  Timers &timers);

void compute_feasible_and_termination(
  const Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  SDP_Solver_Terminate_Reason &terminate_reason, bool &terminate_now);

void compute_dual_residues_and_error(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &y,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
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

SDP_Solver_Terminate_Reason SDP_Solver::run(
  const Environment &env, const Solver_Parameters &parameters,
  const Verbosity &verbosity,
  const boost::property_tree::ptree &parameter_properties,
  const Block_Info &block_info, const SDP &sdp, const El::Grid &grid,
  const std::chrono::time_point<std::chrono::high_resolution_clock> &start_time,
  const fs::path &iterations_json_path, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms)
{
  SDP_Solver_Terminate_Reason terminate_reason(
    SDP_Solver_Terminate_Reason::MaxIterationsExceeded);
  Scoped_Timer solver_timer(timers, "run");
  Scoped_Timer initialize_timer(timers, "initialize");
  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      El::Output(boost::posix_time::second_clock::local_time(),
                 " Start solver iterations");
    }

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

  std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    A_X_inv;

  // BilinearPairingsY is analogous to BilinearPairingsXInv, with
  // X^{-1} -> Y.

  std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    A_Y;
  print_header(verbosity);

  auto psd_sizes(block_info.psd_matrix_block_sizes());
  std::size_t total_psd_rows(
    std::accumulate(psd_sizes.begin(), psd_sizes.end(), size_t(0)));

  Scoped_Timer initialize_bigint_syrk_context_timer(timers,
                                                    "bigint_syrk_context");
  auto bigint_syrk_context = initialize_bigint_syrk_context(
    env, block_info, sdp, parameters.max_shared_memory_bytes,
    verbosity >= Verbosity::debug);
  initialize_bigint_syrk_context_timer.stop();

  initialize_timer.stop();
  auto last_checkpoint_time(std::chrono::high_resolution_clock::now());

  if(El::mpi::Rank() == 0)
    {
      // Copy old iterations.json e.g. to iterations.0.json
      if(fs::exists(iterations_json_path))
        {
          const auto parent_dir = iterations_json_path.parent_path();
          for(size_t index = 0; index < std::numeric_limits<size_t>::max();
              ++index)
            {
              const fs::path backup_path
                = parent_dir
                  / ("iterations." + std::to_string(index) + ".json");
              if(!fs::exists(backup_path))
                {
                  if(verbosity >= Verbosity::debug)
                    El::Output("Move old ", iterations_json_path, " to ",
                               backup_path);
                  fs::rename(iterations_json_path, backup_path);
                  break;
                }
            }
        }

      // Open JSON array
      std::ofstream iterations_json;
      iterations_json.open(iterations_json_path);
      if(iterations_json.good())
        {
          iterations_json << "[";
        }
      else
        {
          if(!iterations_json_path.empty())
            PRINT_WARNING("Cannot write to ", iterations_json_path);
        }
    }

  for(size_t iteration = 1;; ++iteration)
    {
      Scoped_Timer iteration_timer(timers,
                                   "iter_" + std::to_string(iteration));
      if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
        {
          El::Output("Start iteration ", iteration, " at ",
                     boost::posix_time::second_clock::local_time());
        }

      // Prepare graceful exit if any process has received SIGTERM
      {
        El::byte sigterm = env.sigterm_received();
        sigterm = El::mpi::AllReduce(sigterm, El::mpi::LOGICAL_OR,
                                     El::mpi::COMM_WORLD);
        if(sigterm)
          {
            if(El::mpi::Rank() == 0 && !iterations_json_path.empty())
              {
                std::ofstream iterations_json;
                iterations_json.open(iterations_json_path, std::ios::app);
                if(iterations_json.good())
                  iterations_json << "\n]";
              }
            return SDP_Solver_Terminate_Reason::SIGTERM_Received;
          }
      }

      El::byte checkpoint_now(
        std::chrono::duration_cast<std::chrono::seconds>(
          std::chrono::high_resolution_clock::now() - last_checkpoint_time)
          .count()
        >= parameters.checkpoint_interval);
      // Time varies between cores, so follow the decision of the root.
      El::mpi::Broadcast(checkpoint_now, 0, El::mpi::COMM_WORLD);
      if(checkpoint_now == true)
        {
          Scoped_Timer save_timer(timers, "save_checkpoint");
          save_checkpoint(parameters.checkpoint_out, verbosity,
                          parameter_properties);
          last_checkpoint_time = std::chrono::high_resolution_clock::now();
        }
      compute_objectives(sdp, x, y, primal_objective, dual_objective,
                         duality_gap, timers);

      {
        Scoped_Timer cholesky_decomposition_timer(timers,
                                                  "choleskyDecomposition");
        cholesky_decomposition(X, X_cholesky, block_info, "X");
        cholesky_decomposition(Y, Y_cholesky, block_info, "Y");
      }

      compute_bilinear_pairings(block_info, X_cholesky, Y, sdp.bases_blocks,
                                A_X_inv, A_Y, timers);

      compute_dual_residues_and_error(block_info, sdp, y, A_Y, dual_residues,
                                      dual_error, timers);
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
        primal_step_length, dual_step_length, iteration, start_time,
        is_primal_and_dual_feasible, terminate_reason, terminate_now);
      if(terminate_now)
        {
          break;
        }

      El::BigFloat mu, beta_corrector;
      step(parameters, total_psd_rows, is_primal_and_dual_feasible, block_info,
           sdp, grid, X_cholesky, Y_cholesky, A_X_inv, A_Y, primal_residue_p,
           bigint_syrk_context, mu, beta_corrector, primal_step_length,
           dual_step_length, terminate_now, timers, block_timings_ms);

      if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
        {
          El::Print(block_timings_ms, "block_timings, ms:");
          El::Output();
        }
      if(iteration == 1)
        {
          // One the first iteration, matrices may have many zeros, thus
          // block timings may be quite different from the next iterations.
          // Thus, we never want to write first iteration to ck/block_timings.
          if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
            {
              El::Output("block_timings from the first iteration will be "
                         "ignored and removed.");
            }
          block_timings_ms.Empty(false);
        }

      if(terminate_now)
        {
          terminate_reason
            = SDP_Solver_Terminate_Reason::MaxComplementarityExceeded;
          break;
        }
      Scoped_Timer print_iteration_timer(timers, "print_iteration");
      print_iteration(iterations_json_path, iteration, mu, primal_step_length,
                      dual_step_length, beta_corrector, *this,
                      solver_timer.start_time(), verbosity);
    }
  if(El::mpi::Rank() == 0 && !iterations_json_path.empty())
    {
      std::ofstream iterations_json;
      iterations_json.open(iterations_json_path, std::ios::app);
      if(iterations_json.good())
        iterations_json << "\n]";
    }
  return terminate_reason;
}
