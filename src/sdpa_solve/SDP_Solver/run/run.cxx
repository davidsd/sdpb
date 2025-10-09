#include "sdpa_solve/memory_estimates.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "initialize_bigint_syrk_context.hxx"
#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <boost/date_time/posix_time/posix_time.hpp>

// The main solver loop

namespace fs = std::filesystem;

void compute_feasible_and_termination(
  const Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  SDP_Solver_Terminate_Reason &terminate_reason, bool &terminate_now);

namespace Sdpb::Sdpa
{
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
    const El::BigFloat &Q_cond_number,
    const El::BigFloat &max_block_cond_number,
    const std::string &max_block_cond_number_name, const Verbosity &verbosity);

  void compute_objectives(const SDP &sdp, const Primal_Dist_Vector &x,
                          const Block_Diagonal_Matrix &Y,
                          El::BigFloat &primal_objective,
                          El::BigFloat &dual_objective,
                          El::BigFloat &duality_gap, Timers &timers);

  void
  compute_dual_residues_and_error(const SDP &sdp,
                                  const Block_Diagonal_Matrix &Y,
                                  Primal_Dist_Vector &dual_residues,
                                  El::BigFloat &dual_error, Timers &timers);

  void compute_primal_residues_and_error(
    const SDP &sdp, const Primal_Dist_Vector &x,
    const Block_Diagonal_Matrix &X, Block_Diagonal_Matrix &primal_residues,
    El::BigFloat &primal_error, Timers &timers);

  namespace
  {
    // Memory allocated in initialize_P(), see in step/compute_S.cxx
    // Includes:
    // - Size of matrix P (Block_Matrix part on current node)
    // - Memory usage for processing the largest block in a group:
    //   - F_p L_Y block
    //   - G_p = L_X_inv F_p L_Y.
    //   - Extra memory usage inside El::Trsm() call.
    // NB: returns all memory for MPI group from group_rank = 0
    // Other ranks return 0, so that mpi::Reduce() returns correct total memory.
    size_t
    get_initialize_P_bytes(const Block_Info &block_info, const El::Grid &grid)
    {
      const auto &group_comm = block_info.mpi_comm.value;
      const size_t primal_dimension = block_info.primal_dimension;

      if(group_comm.Rank() != 0)
        return 0;

      // Height of P submatrix stored on current MPI group
      size_t P_group_height = 0;
      const size_t P_width = primal_dimension;

      // Choose the largest block in the group
      size_t max_block_dim = 0;
      for(const auto &block_index : block_info.block_indices)
        {
          const size_t dim = block_info.block_dimensions.at(block_index);
          P_group_height += dim * dim;
          max_block_dim = std::max(max_block_dim, dim);
        }

      const size_t P_size = P_group_height * P_width;
      const size_t P_bytes = bigfloat_bytes() * P_size;

      const size_t FY_block_vertical_bytes
        = bigfloat_bytes() * max_block_dim * max_block_dim * primal_dimension;
      const size_t G_block_horizontal_bytes = FY_block_vertical_bytes;

      const size_t trsm_bytes
        = get_trsm_bytes(max_block_dim, max_block_dim * primal_dimension,
                         grid.Height(), grid.Width());
      return P_bytes + FY_block_vertical_bytes + G_block_horizontal_bytes
             + trsm_bytes;
    }

    // Estimate how many BigFloats will be allocated by SDPB on the current node,
    // (including what's already allocated, e.g. SDP)
    size_t get_required_nonshared_memory_per_node_bytes(
      const Environment &env, const Block_Info &block_info, const SDP &sdp,
      const El::Grid &grid, const SDP_Solver &solver,
      const Verbosity verbosity)
    {
      const auto &node_comm = env.comm_shared_mem;

      const auto node_reduce = [&node_comm](const size_t value) -> size_t {
        return El::mpi::Reduce(value, 0, node_comm);
      };

      // Block_Diagonal_Matrix X, Y, primal_residues, X_chol, Y_chol, dX, dY, R, Z
      const size_t X_size = node_reduce(get_matrix_size_local(solver.X));
      const size_t X_bytes = node_reduce(get_allocated_bytes(solver.X));

      // Primal_Dist_Vector x
      const size_t x_bytes = node_reduce(get_allocated_bytes(solver.x));

      // DistMatrix S ~ MxM, distributed over all nodes.
      const size_t S_size = node_reduce(get_S_size_local(sdp));
      const size_t S_bytes = node_reduce(get_S_allocated_bytes(sdp));

      // SDP struct
      // const size_t SDP_size = node_reduce(get_SDP_size_local(sdp));
      const size_t SDP_bytes = node_reduce(get_allocated_bytes(sdp));

      // SDP_Solver members
      const size_t SDP_solver_bytes = node_reduce(get_allocated_bytes(solver));

      const size_t initialize_P_bytes
        = node_reduce(get_initialize_P_bytes(block_info, grid));

      // run(): X_cholesky, Y_cholesky
      size_t SDP_Solver_run_bytes = 2 * X_bytes;
      // step(): dx, dX, dY, S
      SDP_Solver_run_bytes += x_bytes + 2 * X_bytes + S_bytes;
      // step():
      SDP_Solver_run_bytes += std::max(
        // Temporary allocations inside compute_S()
        initialize_P_bytes,
        // Allocations after compute_S():
        // step(): minus_XY
        // compute_search_direction(): R, Z
        3 * X_bytes);

      // We will use only result on rank=0
      if(node_comm.Rank() != 0)
        return 0;

      // Calculate mem_required_size

      // initial_node_mem_used() is RAM allocated at SDPB start.
      // This could be important: e.g. on 128 cores (Expanse HPC) it is ~26GB
      const size_t mem_required_bytes = env.initial_node_mem_used() + SDP_bytes
                                        + SDP_solver_bytes
                                        + SDP_Solver_run_bytes;

      if(verbosity >= Verbosity::debug)
        {
          std::ostringstream ss;
          El::BuildStream(
            ss, "node=", env.node_index(),
            " matrix sizes and memory estimates: ", "\n\t#(X) = ", X_size,
            "\n\t#(S) = ", S_size, "\n memory estimates: ",
            "\n\tBigFloat size: ", pretty_print_bytes(bigfloat_bytes(), true),
            "\n\t#(SDP) = ", pretty_print_bytes(SDP_bytes, true),
            "\n\tInitial MemUsed (at SDPB start) = ",
            pretty_print_bytes(env.initial_node_mem_used(), true),
            "\n\tTotal non-shared memory estimate: ",
            pretty_print_bytes(mem_required_bytes, true));
          El::Output(ss.str());
        }

      return mem_required_bytes;
    }

    size_t
    get_max_shared_memory_bytes(const size_t default_max_shared_memory_bytes,
                                const Environment &env,
                                const Block_Info &block_info, const SDP &sdp,
                                const El::Grid &grid, const SDP_Solver &solver,
                                const Verbosity verbosity)
    {
      // If user sets --maxSharedMemory limit manually, we use it.
      // Otherwise, we calculate the limit automatically.
      if(default_max_shared_memory_bytes != 0)
        return default_max_shared_memory_bytes;
      const size_t nonshared_memory_required_per_node_bytes
        = get_required_nonshared_memory_per_node_bytes(
          env, block_info, sdp, grid, solver, verbosity);
      return get_max_shared_memory_bytes(
        nonshared_memory_required_per_node_bytes, env, verbosity);
    }
  }

  SDP_Solver_Terminate_Reason SDP_Solver::run(
    const Environment &env, const Solver_Parameters &parameters,
    const Verbosity &verbosity,
    const boost::property_tree::ptree &parameter_properties,
    const Block_Info &block_info, const SDP &sdp, const El::Grid &grid,
    const std::chrono::time_point<std::chrono::high_resolution_clock>
      &start_time,
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
    if(verbosity >= Verbosity::debug)
      {
        print_allocation_message_per_node(env, "X_cholesky",
                                          get_allocated_bytes(X_cholesky));
        print_allocation_message_per_node(env, "Y_cholesky",
                                          get_allocated_bytes(Y_cholesky));
      }

    size_t total_psd_rows
      = std::accumulate(block_info.block_dimensions.begin(),
                        block_info.block_dimensions.end(), size_t(0));

    Scoped_Timer initialize_bigint_syrk_context_timer(timers,
                                                      "bigint_syrk_context");

    auto max_shared_memory_bytes
      = get_max_shared_memory_bytes(parameters.max_shared_memory_bytes, env,
                                    block_info, sdp, grid, *this, verbosity);
    auto bigint_syrk_context = initialize_bigint_syrk_context(
      env, block_info, sdp, max_shared_memory_bytes, verbosity);
    initialize_bigint_syrk_context_timer.stop();

    initialize_timer.stop();
    auto last_checkpoint_time(std::chrono::high_resolution_clock::now());

    const bool write_iterations_json_and_c_minus_By
      = !iterations_json_path.empty();
    // TODO: we don't have (c - B.y) anymore
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
                      iterations_json.open(iterations_json_path,
                                           std::ios::app);
                      if(iterations_json.good())
                        iterations_json << "\n]";
                    }
                  // auto c_minus_By_path
                  //   = c_minus_By_dir_path / ("c_minus_By.json");
                  // save_c_minus_By(c_minus_By_path, block_info, sdp, y, verbosity,
                  //                 timers);
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
                // auto c_minus_By_path
                //   = c_minus_By_dir_path
                //     / ("c_minus_By." + std::to_string(iteration) + ".json");
                // save_c_minus_By(c_minus_By_path, block_info, sdp, y, verbosity,
                //                 timers);
              }
          }
        compute_objectives(sdp, x, Y, primal_objective, dual_objective,
                           duality_gap, timers);

        {
          Scoped_Timer cholesky_decomposition_timer(timers,
                                                    "choleskyDecomposition");
          cholesky_decomposition(X, X_cholesky, block_info.block_indices, "X");
          cholesky_decomposition(Y, Y_cholesky, block_info.block_indices, "Y");
        }

        compute_dual_residues_and_error(sdp, Y, dual_residues, dual_error,
                                        timers);
        compute_primal_residues_and_error(sdp, x, X, primal_residues,
                                          primal_error, timers);

        bool terminate_now, is_primal_and_dual_feasible;
        compute_feasible_and_termination(
          parameters, primal_error, dual_error, duality_gap,
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
             Y_cholesky, bigint_syrk_context, mu, beta_corrector,
             primal_step_length, dual_step_length, terminate_now, timers,
             block_timings_ms, Q_cond_number, max_block_cond_number,
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
        print_iteration(
          iterations_json_path, iteration, mu, primal_step_length,
          dual_step_length, beta_corrector, *this, solver_timer.start_time(),
          iteration_timer.start_time(), Q_cond_number, max_block_cond_number,
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
      }
    return terminate_reason;
  }
}
