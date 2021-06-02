#pragma once
// Parameters to control the behavior of an SDPSolver. See the manual
// for a detailed description of each.
//

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>

#include <iostream>

struct Solver_Parameters
{
  int64_t max_iterations, max_runtime, checkpoint_interval;
  bool find_primal_feasible, find_dual_feasible, detect_primal_feasible_jump,
    detect_dual_feasible_jump;
  size_t precision;

  El::BigFloat duality_gap_threshold, primal_error_threshold,
    dual_error_threshold, initial_matrix_scale_primal,
    initial_matrix_scale_dual, feasible_centering_parameter,
    infeasible_centering_parameter, step_length_reduction, max_complementarity,
    min_primal_step, min_dual_step;

  boost::filesystem::path checkpoint_in, checkpoint_out;
  Solver_Parameters() = default;
  boost::program_options::options_description options();
};

std::ostream &operator<<(std::ostream &os, const Solver_Parameters &p);

boost::property_tree::ptree to_property_tree(const Solver_Parameters &p);
