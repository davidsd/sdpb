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
#include <vector>
struct Dynamical_Solver_Parameters
{
  Solver_Parameters solver_parameters;  

  El::BigFloat alpha,
               //duality_gap_threshold, primal_error_threshold, dual_error_threshold ,
               //primal_obj_threshold, dual_obj_threshold,
               update_sdp_threshold_max, update_sdp_threshold_min,
               find_boundary_obj_threshold;

  El::BigFloat updateSDP_dualityGapThreshold, dualityGap_upper_limit;

  int n_external_parameters;
  size_t total_iterations;
  boost::filesystem::path new_sdp_path;

  boost::filesystem::path old_schur_path; 
  boost::filesystem::path old_sdp_path;
  bool external_corrector_Q;

  std::vector<El::BigFloat> bounding_box_max;
  std::vector<El::BigFloat> bounding_box_min;

  bool find_boundary, fix_ext_param_direction;
  std::vector<El::BigFloat> search_direction; 

  bool use_exact_hessian; 
  std::vector<El::BigFloat> prev_grad;
  std::vector<El::BigFloat> prev_step;
  std::vector<El::BigFloat> hess_BFGS;

  bool printMore;
  bool use_Hmixed_for_BFGS;
  El::BigFloat centeringRThreshold;


  // new strategy parameters
  El::BigFloat beta_scan_min;
  El::BigFloat beta_scan_max;
  El::BigFloat beta_scan_step;

  El::BigFloat step_min_threshold;
  El::BigFloat step_max_threshold;
  El::BigFloat lagrangian_muI_shift;

  El::BigFloat finite_dGap_target;

  int max_climbing;
  El::BigFloat beta_climbing;
  El::BigFloat beta_for_mu_logdetX;
  El::BigFloat BFGS_partial_update_reduction;

  bool gradientWithLogDetX, stickToGCPQ, optimalbetaQ;

  El::BigFloat navigatorValueShift;
  bool navigatorAutomaticShiftQ;

  bool climbingRecomputeExtParamQ;

  bool returnCheckpointOnLCP;

  Dynamical_Solver_Parameters() = default;
  boost::program_options::options_description options();

};

std::ostream &operator<<(std::ostream &os, const Dynamical_Solver_Parameters &p);

boost::property_tree::ptree to_property_tree(const Dynamical_Solver_Parameters &p);


struct PrecParameters
{
	size_t prec;
	PrecParameters(int argc, char *argv[]);
};


