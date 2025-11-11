#pragma once

#include "SDP.hxx"
#include "SDP_Solver/run/step/Compute_S_Context.hxx"
#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/SDP_Solver_Terminate_Reason.hxx"

#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/Solver_Parameters.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>

namespace Sdpb::Sdpa
{
  // SDPSolver contains the data structures needed during the running of
  // the interior point algorithm.  Each structure is allocated when an
  // SDPSolver is initialized, and reused in each iteration.
  //
  struct SDP_Solver
  {
    using Block_Info_Type = Block_Info;
    using SDP_Type = SDP;

    // a Vector of length m
    // Duplicated among all blocks.
    // TODO: why the type is different from the vector SDP::primal_objective_c?
    Primal_Dist_Vector x;

    // a Block_Diagonal_Matrix with block sizes given by
    // sdp.psdMatrixBlockDims()
    Block_Diagonal_Matrix X;

    // a Block_Diagonal_Matrix with the same structure as X
    Block_Diagonal_Matrix Y;

    /********************************************/
    // Solver status

    El::BigFloat primal_objective, // c . x
      dual_objective,              // Tr (F_0 Y)
      duality_gap;                 // normalized difference of objectives

    // Discrepancy in the primal equality constraints, a
    // Block_Diagonal_Matrix with the same structure as X, called 'P' in
    // the manual:
    //
    //   PrimalResidues = \sum_i F_i x_i - F_0 - X
    //
    Block_Diagonal_Matrix primal_residues;

    // |P|
    El::BigFloat primal_error;

    // Discrepancy in the dual equality constraints, a Vector of length
    // P, called 'd' in the manual:
    //
    //   dualResidues = c - Tr(F_* Y)
    //
    Primal_Dist_Vector dual_residues;
    El::BigFloat dual_error; // maxAbs(dualResidues)
    El::BigFloat R_error;    // maxAbs(R = mu * I - XY)

    int64_t current_generation;
    boost::optional<int64_t> backup_generation;

    SDP_Solver(const Solver_Parameters &parameters, const Verbosity &verbosity,
               const bool &require_initial_checkpoint,
               const Block_Info &block_info, const El::Grid &grid);

    SDP_Solver_Terminate_Reason
    run(const Environment &env, const Solver_Parameters &parameters,
        const Verbosity &verbosity,
        const boost::property_tree::ptree &parameter_properties,
        const Block_Info &block_info, const SDP &sdp, const El::Grid &grid,
        const std::chrono::time_point<std::chrono::high_resolution_clock>
          &start_time,
        const std::filesystem::path &iterations_json_path, Timers &timers,
        El::Matrix<int32_t> &block_timings_ms);

    void step(const Environment &env, const Solver_Parameters &parameters,
              const Verbosity &verbosity, const std::size_t &total_psd_rows,
              const bool &is_primal_and_dual_feasible,
              const Block_Info &block_info, const SDP &sdp,
              const El::Grid &grid, const Block_Diagonal_Matrix &X_cholesky,
              const Block_Diagonal_Matrix &Y_cholesky,
              Compute_S_Context &compute_S_context, El::BigFloat &mu,
              El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
              El::BigFloat &dual_step_length, bool &terminate_now,
              Timers &timers, El::Matrix<int32_t> &block_timings_ms,
              El::BigFloat &Q_cond_number, El::BigFloat &max_block_cond_number,
              std::string &max_block_cond_number_name);

    void
    save_checkpoint(const std::filesystem::path &checkpoint_directory,
                    const Verbosity &verbosity,
                    const boost::property_tree::ptree &parameter_properties);
    bool
    load_checkpoint(const std::filesystem::path &checkpoint_directory,
                    const Block_Info &block_info, const Verbosity &verbosity,
                    const bool &require_initial_checkpoint);
  };
}
