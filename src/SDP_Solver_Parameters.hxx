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

  El::BigFloat duality_gap_threshold_elemental, primal_error_threshold_elemental,
    dual_error_threshold_elemental, initial_matrix_scale_primal_elemental,
    initial_matrix_scale_dual_elemental,
    feasible_centering_parameter_elemental,
    infeasible_centering_parameter_elemental, step_length_reduction_elemental,
    max_complementarity_elemental;

  friend std::ostream &
  operator<<(std::ostream &os, const SDP_Solver_Parameters &p);
};
