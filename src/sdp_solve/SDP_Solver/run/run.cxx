#include "save_c_minus_By.hxx"
#include "sdp_solve/memory_estimates.hxx"
#include "bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "get_syrk_Q_config.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <boost/date_time/posix_time/posix_time.hpp>

// The main solver loop

namespace fs = std::filesystem;

void print_header(const Verbosity &verbosity);
void print_iteration(
  const fs::path &iterations_json_path, const int &iteration,
  const El::BigFloat &mu, const El::BigFloat &primal_step_length,
  const El::BigFloat &dual_step_length, const El::BigFloat &beta_corrector,
  const SDP_Solver &sdp_solver,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &iteration_start_time,
  const El::BigFloat &Q_cond_number, const El::BigFloat &max_block_cond_number,
  const std::string &max_block_cond_number_name, const Verbosity &verbosity);

void compute_objectives(const SDP &sdp, const Block_Vector &x,
                        const Block_Vector &y, El::BigFloat &primal_objective,
                        El::BigFloat &dual_objective,
                        El::BigFloat &duality_gap, Timers &timers);

void compute_bilinear_pairings(
  const Block_Info &block_info, const Paired_Block_Diagonal_Matrix &X_cholesky,
  const Paired_Block_Diagonal_Matrix &Y,
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
  const Paired_Block_Diagonal_Matrix &X,
  Paired_Block_Diagonal_Matrix &primal_residues, El::BigFloat &primal_error_P,
  Timers &timers);

void compute_primal_residues_and_error_p_b_Bx(const Block_Info &block_info,
                                              const SDP &sdp,
                                              const Block_Vector &x,
                                              Block_Vector &primal_residue_p,
                                              El::BigFloat &primal_error_p);

namespace
{
  BigInt_Shared_Memory_Syrk_Context create_syrk_Q_context_and_print_memory(
    const Environment &env, const Block_Info &block_info, const El::Grid &grid,
    const SDP &sdp, const SDP_Solver &solver,
    const Solver_Parameters &parameters, const Verbosity verbosity)
  {
    // Estimate how many BigFloats will be allocated by SDPB on the current node,
    // (including what's already allocated, e.g. SDP)
    const auto &node_comm = env.comm_shared_mem;

    const auto node_reduce = [&node_comm](const size_t value) -> size_t {
      return El::mpi::AllReduce(value, node_comm);
    };

    // X, Y, primal_residues, X_chol, Y_chol, dX, dY, R, Z
    const size_t X_size = node_reduce(get_matrix_size_local(solver.X));
    const size_t X_bytes = node_reduce(get_allocated_bytes(solver.X));

    // Bilinear pairing blocks - A_X_inv, A_Y
    const size_t A_X_inv_size
      = node_reduce(get_A_X_size_local(block_info, sdp));
    const size_t A_X_inv_bytes = A_X_inv_size * bigfloat_bytes();

    // schur_complement, schur_complement_cholesky
    const size_t schur_complement_size
      = node_reduce(get_schur_complement_size_local(block_info));
    const size_t schur_complement_bytes
      = schur_complement_size * bigfloat_bytes();
    const size_t schur_complement_cholesky_bytes = schur_complement_bytes;

    // Memory overhead for computing L^{-1}
    // (Trsm calls in initialize_schur_off_diagonal).
    // Sum over MPI groups the Trsm memory for the largest block of B in group.
    const size_t initialize_schur_off_diagonal_trsm_bytes = node_reduce([&] {
      size_t trsm_bytes = 0;
      if(block_info.mpi_comm.value.Rank() == 0)
        {
          for(auto &m : sdp.free_var_matrix.blocks)
            {
              trsm_bytes = std::max(
                trsm_bytes, get_trsm_bytes(m.Height(), m.Width(),
                                           grid.Height(), grid.Width()));
            }
        }
      return trsm_bytes;
    }());

    // #(B) = PxN
    // sdp.free_var_matrix, schur_off_diagonal
    const size_t B_size = node_reduce(get_B_size_local(sdp));
    const size_t B_bytes
      = node_reduce(get_allocated_bytes(sdp.free_var_matrix));

    // #Q = NxN, distributed over all nodes.
    const size_t Q_size = node_reduce(get_Q_size_local(sdp));
    const size_t Q_bytes = Q_size * bigfloat_bytes();

    // SDP struct
    const size_t SDP_bytes = node_reduce(get_allocated_bytes(sdp));

    // SDP_Solver struct
    const size_t SDP_Solver_bytes = node_reduce(get_allocated_bytes(solver));

    // SDP_Solver::run()
    size_t SDP_Solver_run_bytes = 0;

    const auto mem_required_bytes
      = [&] { return SDP_bytes + SDP_Solver_bytes + SDP_Solver_run_bytes; };

    // X_cholesky, Y_cholesky, primal_residues, dX, dY
    SDP_Solver_run_bytes += 5 * X_bytes;

    // A_X_inv and A_Y
    SDP_Solver_run_bytes += 2 * A_X_inv_bytes;

    // schur_complement_cholesky
    SDP_Solver_run_bytes += schur_complement_cholesky_bytes;

    // schur_off_diagonal = L^{-1} B
    SDP_Solver_run_bytes += B_bytes;
    // Q = NxN
    SDP_Solver_run_bytes += Q_bytes;

    const auto cfg = get_syrk_Q_config(
      block_info, sdp, parameters.max_memory, parameters.max_shared_memory,
      mem_required_bytes() + std::max(schur_complement_bytes, 3 * X_bytes));

    // Add schur_complement + (syrk_S() or Trsm() temporary allocations)
    // from initialize_schur_complement_solver()
    // or XY,R,Z from compute_search_direction().
    // (when allocations do not coexist, we choose maximum size instead of adding both)
    // TODO implement RAII-like simulation of memory allocations/deallocations,
    // to automatically track and print max memory usage.
    const size_t initialize_schur_complement_solver_local_bytes
      = schur_complement_bytes
        + std::max(cfg.node_local_bytes(),
                   initialize_schur_off_diagonal_trsm_bytes);
    const size_t initialize_schur_complement_solver_total_bytes
      = initialize_schur_complement_solver_local_bytes
        + cfg.node_shmem_bytes();
    SDP_Solver_run_bytes
      += std::max(initialize_schur_complement_solver_local_bytes, 3 * X_bytes);
    // Shared memory windows are allocated in the beginning, so they should be always included.
    SDP_Solver_run_bytes += cfg.node_shmem_bytes();

    if(node_comm.Rank() == 0 && verbosity >= Verbosity::debug)
      {
        // Print memory estimates

        std::vector<std::pair<std::string, size_t>> num_elements_per_category{
          {"X", X_size},
          {"A_X_inv", A_X_inv_size},
          {"schur_complement", schur_complement_size},
          {"B", B_size},
          {"Q", Q_size},
        };

        std::vector<std::pair<std::string, size_t>> bytes_per_category{
          {"Initial MemAvailable (at SDPB start)",
           env.initial_node_mem_available()},
          {"BigFloat size", bigfloat_bytes()},
          {"Total SDPB memory estimate", mem_required_bytes()},
          {"Shared memory estimate", cfg.node_shmem_bytes()},
          {"\tSDP", SDP_bytes},
          {"\tSDP_Solver", SDP_Solver_bytes},
          {"\tSDP_Solver::run()", SDP_Solver_run_bytes},
          {"\t\tinitialize_schur_complement_solver()",
           initialize_schur_complement_solver_total_bytes},
          {"\t\t\tschur_complement", schur_complement_bytes},
          {"\t\t\tEl::Trsm()", initialize_schur_off_diagonal_trsm_bytes},
          {"\t\t\tQ", Q_bytes},
          {"\t\t\tshared memory", cfg.node_shmem_bytes()},
        };

        std::ostringstream ss;
        El::BuildStream(ss, "node=", env.node_index(),
                        " matrix sizes and memory estimates: ");

        for(const auto &[name, size] : num_elements_per_category)
          {
            El::BuildStream(ss, "\n\t#(", name, ") = ", size, " elements");
          }
        for(const auto &[name, bytes] : bytes_per_category)
          {
            El::BuildStream(ss, "\n\t", name, ": ",
                            pretty_print_bytes(bytes, true));
          }
        El::BuildStream(
          ss, "\n\tShared memory configuration: ",
          "\n\t\tsyrk_Q() input split factor: ", cfg.input_split_factor,
          "\n\t\tsyrk_Q() output split factor: ", cfg.output_split_factor);
        El::Output(ss.str());
      }

    return BigInt_Shared_Memory_Syrk_Context(
      cfg, block_info.node_group_index(), block_info.block_indices, verbosity);
  }
}

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

  auto X_cholesky(X), Y_cholesky(X);
  if(verbosity >= Verbosity::debug)
    {
      print_allocation_message_per_node(env, "X_cholesky",
                                        get_allocated_bytes(X_cholesky));
      print_allocation_message_per_node(env, "Y_cholesky",
                                        get_allocated_bytes(Y_cholesky));
    }

  // Bilinear pairings needed for computing the Schur complement
  // matrix.  For example,
  //
  //   BilinearPairingsXInv.blocks[b].elt(
  //     (d_j+1) s + k1,
  //     (d_j+1) r + k2
  //   ) = v_{b,k1}^T (X.blocks[b]^{-1})^{(s,r)} v_{b,k2}
  //
  //     0 <= k1,k2 < block_info.num_points[j] = d_j + 1
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

  auto psd_sizes(block_info.psd_matrix_block_sizes());
  const size_t total_psd_rows
    = std::accumulate(psd_sizes.begin(), psd_sizes.end(), 0);

  Scoped_Timer initialize_bigint_syrk_context_timer(timers,
                                                    "bigint_syrk_context");
  auto bigint_syrk_context = create_syrk_Q_context_and_print_memory(
    env, block_info, grid, sdp, *this, parameters, verbosity);
  initialize_bigint_syrk_context_timer.stop();

  initialize_timer.stop();
  auto last_checkpoint_time(std::chrono::high_resolution_clock::now());

  const bool write_iterations_json_and_c_minus_By
    = !iterations_json_path.empty();
  const auto c_minus_By_dir_path
    = write_iterations_json_and_c_minus_By
        ? iterations_json_path.parent_path() / "c_minus_By"
        : fs::path();
  if(El::mpi::Rank() == 0 && write_iterations_json_and_c_minus_By)
    {
      // Copy old out/iterations.json e.g. to out/iterations.0.json,
      // Do the same for out/c_minus_By/
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

                  // move old c_minus_By
                  if(fs::exists(c_minus_By_dir_path))
                    {
                      auto c_minus_By_dir_backup_path = c_minus_By_dir_path;
                      c_minus_By_dir_backup_path.replace_extension(
                        "." + std::to_string(index));
                      if(verbosity >= Verbosity::debug)
                        El::Output("Move old ", c_minus_By_dir_path, " to ",
                                   c_minus_By_dir_backup_path);
                      if(fs::exists(c_minus_By_dir_backup_path))
                        PRINT_WARNING(c_minus_By_dir_backup_path,
                                      " exists and will be overwritten");
                      fs::rename(c_minus_By_dir_path,
                                 c_minus_By_dir_backup_path);
                    }

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

  print_header(verbosity);
  for(size_t iteration = 1;; ++iteration)
    {
      Scoped_Timer iteration_timer(timers,
                                   "iter_" + std::to_string(iteration));
      if(verbosity >= Verbosity::trace && El::mpi::Rank() == 0)
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
            if(write_iterations_json_and_c_minus_By)
              {
                if(El::mpi::Rank() == 0)
                  {
                    std::ofstream iterations_json;
                    iterations_json.open(iterations_json_path, std::ios::app);
                    if(iterations_json.good())
                      iterations_json << "\n]";
                  }
                auto c_minus_By_path
                  = c_minus_By_dir_path / ("c_minus_By.json");
                save_c_minus_By(c_minus_By_path, block_info, sdp, y, verbosity,
                                timers);
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

          if(write_iterations_json_and_c_minus_By)
            {
              auto c_minus_By_path
                = c_minus_By_dir_path
                  / ("c_minus_By." + std::to_string(iteration) + ".json");
              save_c_minus_By(c_minus_By_path, block_info, sdp, y, verbosity,
                              timers);
            }
        }
      compute_objectives(sdp, x, y, primal_objective, dual_objective,
                         duality_gap, timers);

      {
        Scoped_Timer cholesky_decomposition_timer(timers,
                                                  "choleskyDecomposition");
        cholesky_decomposition(X, X_cholesky, block_info.block_indices, "X");
        cholesky_decomposition(Y, Y_cholesky, block_info.block_indices, "Y");
      }

      compute_bilinear_pairings(block_info, X_cholesky, Y, sdp.bases_blocks,
                                A_X_inv, A_Y, timers);
      if(verbosity >= Verbosity::trace)
        {
          print_allocation_message_per_node(env, "A_X_inv",
                                            get_allocated_bytes(A_X_inv));
          print_allocation_message_per_node(env, "A_Y",
                                            get_allocated_bytes(A_Y));
        }

      compute_dual_residues_and_error(block_info, sdp, y, A_Y, dual_residues,
                                      dual_error, timers);
      compute_primal_residues_and_error_P_Ax_X(
        block_info, sdp, x, X, primal_residues, primal_error_P, timers);

      // use y to set the sizes of primal_residue_p.  The data is
      // overwritten in compute_primal_residues_and_error_p.
      Block_Vector primal_residue_p(y);
      compute_primal_residues_and_error_p_b_Bx(
        block_info, sdp, x, primal_residue_p, primal_error_p);
      if(verbosity >= Verbosity::trace)
        {
          print_allocation_message_per_node(
            env, "primal_residue_p", get_allocated_bytes(primal_residue_p));
        }

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
      El::BigFloat Q_cond_number;
      El::BigFloat max_block_cond_number;
      std::string max_block_cond_number_name;
      step(env, parameters, verbosity, total_psd_rows,
           is_primal_and_dual_feasible, block_info, sdp, grid, X_cholesky,
           Y_cholesky, A_X_inv, A_Y, primal_residue_p, bigint_syrk_context, mu,
           beta_corrector, primal_step_length, dual_step_length, terminate_now,
           timers, block_timings_ms, Q_cond_number, max_block_cond_number,
           max_block_cond_number_name);

      if(verbosity >= Verbosity::trace && El::mpi::Rank() == 0)
        {
          El::Print(block_timings_ms, "block_timings, ms:");
          El::Output();
        }
      if(iteration == 1)
        {
          // One the first iteration, matrices may have many zeros, thus
          // block timings may be quite different from the next iterations.
          // Thus, we never want to write first iteration to ck/block_timings.
          if(verbosity >= Verbosity::trace && El::mpi::Rank() == 0)
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
                      solver_timer.start_time(), iteration_timer.start_time(),
                      Q_cond_number, max_block_cond_number,
                      max_block_cond_number_name, verbosity);
    }

  if(write_iterations_json_and_c_minus_By)
    {
      if(El::mpi::Rank() == 0)
        {
          std::ofstream iterations_json;
          iterations_json.open(iterations_json_path, std::ios::app);
          if(iterations_json.good())
            iterations_json << "\n]";
        }
      auto c_minus_By_path = c_minus_By_dir_path / "c_minus_By.json";
      save_c_minus_By(c_minus_By_path, block_info, sdp, y, verbosity, timers);
    }
  return terminate_reason;
}
