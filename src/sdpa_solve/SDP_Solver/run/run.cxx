#include "initialize_compute_S_context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpa_solve/memory_estimates.hxx"

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

    Scoped_Timer initialize_compute_S_context_timer(timers,
                                                    "compute_S_context");

    auto cfg = Solver_Run_Config::create(env, block_info, sdp, *this,
                                         parameters, verbosity);
    auto compute_S_context
      = Compute_S_Context::create(block_info, cfg, verbosity);
    initialize_compute_S_context_timer.stop();

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
             Y_cholesky, compute_S_context, mu, beta_corrector,
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
