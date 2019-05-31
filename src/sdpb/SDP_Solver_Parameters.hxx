#pragma once
// Parameters to control the behavior of an SDPSolver. See the manual
// for a detailed description of each.
//

#include "Verbosity.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>

#include <iostream>

class SDP_Solver_Parameters
{
public:
  int64_t max_iterations, max_runtime, checkpoint_interval;
  bool no_final_checkpoint, find_primal_feasible, find_dual_feasible,
    detect_primal_feasible_jump, detect_dual_feasible_jump, write_matrices;
  size_t precision, procs_per_node, proc_granularity;
  Verbosity verbosity;

  El::BigFloat duality_gap_threshold, primal_error_threshold,
    dual_error_threshold, initial_matrix_scale_primal,
    initial_matrix_scale_dual, feasible_centering_parameter,
    infeasible_centering_parameter, step_length_reduction, max_complementarity;

  boost::filesystem::path sdp_directory, out_directory, checkpoint_in,
    checkpoint_out, param_file;

  SDP_Solver_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !sdp_directory.empty(); }
};

std::ostream &operator<<(std::ostream &os, const SDP_Solver_Parameters &p);
