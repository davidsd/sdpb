#pragma once
// Parameters to control the behavior of an SDPSolver. See the manual
// for a detailed description of each.
//

#include <El.hpp>

#include <iostream>

class SDP_Solver_Parameters
{
public:
  int max_iterations;
  int max_runtime;
  int checkpoint_interval;
  bool no_final_checkpoint;
  bool find_primal_feasible;
  bool find_dual_feasible;
  bool detect_primal_feasible_jump;
  bool detect_dual_feasible_jump;
  int precision;
  int procs_per_node;

  El::BigFloat duality_gap_threshold, primal_error_threshold,
    dual_error_threshold, initial_matrix_scale_primal,
    initial_matrix_scale_dual, feasible_centering_parameter,
    infeasible_centering_parameter, step_length_reduction, max_complementarity;

  friend std::ostream &
  operator<<(std::ostream &os, const SDP_Solver_Parameters &p);
};
