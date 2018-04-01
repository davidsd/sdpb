#pragma once
// Parameters to control the behavior of an SDPSolver. See the manual
// for a detailed description of each.
//

#include "types.hxx"

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
  Real duality_gap_threshold;
  Real primal_error_threshold;
  Real dual_error_threshold;
  Real initial_matrix_scale_primal;
  Real initial_matrix_scale_dual;
  Real feasible_centering_parameter;
  Real infeasible_centering_parameter;
  Real step_length_reduction;
  Real max_complementarity;
  // Set the precision of all Real parameters to equal 'precision'.
  // This is necessary because 'precision' might be set (via the
  // command line or a file) after initializing other parameters.
  //
  void resetPrecision()
  {
    duality_gap_threshold.set_prec(precision);
    primal_error_threshold.set_prec(precision);
    dual_error_threshold.set_prec(precision);
    initial_matrix_scale_primal.set_prec(precision);
    initial_matrix_scale_dual.set_prec(precision);
    feasible_centering_parameter.set_prec(precision);
    infeasible_centering_parameter.set_prec(precision);
    step_length_reduction.set_prec(precision);
    max_complementarity.set_prec(precision);
  }

  friend std::ostream &
  operator<<(std::ostream &os, const SDP_Solver_Parameters &p);
};
