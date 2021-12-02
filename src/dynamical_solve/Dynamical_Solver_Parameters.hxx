#pragma once
// Parameters to control the behavior of an SDPSolver. See the manual
// for a detailed description of each.
//

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>

#include "../sdp_solve/Solver_Parameters.hxx"
#include <iostream>

struct Dynamical_Solver_Parameters
{
  Solver_Parameters solver_parameters;  

  El::BigFloat alpha,update_sdp_threshold_max;
               //duality_gap_threshold, primal_error_threshold, dual_error_threshold ,
               //primal_obj_threshold, dual_obj_threshold,
               //update_sdp_threshold_max, update_sdp_threshold_min;

  int n_external_parameters;
  boost::filesystem::path new_sdp_path;

  Dynamical_Solver_Parameters() = default;
  boost::program_options::options_description options();

};

std::ostream &operator<<(std::ostream &os, const Dynamical_Solver_Parameters &p);

boost::property_tree::ptree to_property_tree(const Dynamical_Solver_Parameters &p);
