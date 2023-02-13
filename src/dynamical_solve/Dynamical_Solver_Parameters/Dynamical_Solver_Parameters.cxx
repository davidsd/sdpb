#include "../Dynamical_Solver_Parameters.hxx"

#include <boost/program_options.hpp>
#include <vector>

boost::program_options::options_description Dynamical_Solver_Parameters::options()
{
  boost::program_options::options_description result("Dynamic Solver parameters");
  // We set default parameters using El::BigFloat("1e-10",10)
  // rather than a straight double precision 1e-10 so that results
  // are more reproducible at high precision.  Using double
  // precision defaults results in differences of about 1e-15 in
  // primalObjective after one step.

  result.add_options()("newSdpDirs", 
                              boost::program_options::value<boost::filesystem::path>(&new_sdp_path)->required(),
                              "Directory containing the preprocessed SDP data files around the center SDP in external parameter space.");  
 
  result.add_options()("stepSizeAlpha", 
                              boost::program_options::value<El::BigFloat>(&alpha)->default_value(1),
                              "Step size in the external-parameter space to generate the new SDP data files. "
                              "The default value is set to 1.");
  result.add_options()("numExternalParams",
                              boost::program_options::value<int>(&n_external_parameters)->default_value(0),
                              "The number of external parameters to be varied in each iteration of the dynamical SDP. "
                              "The default value is set to 0.");
   result.add_options()("prevTotalIterations", 
                              boost::program_options::value<size_t>(&total_iterations)->default_value(0),
                              "The number of total iterations finished before this solver starts to run. "
                              "The default value is set to 0.");
  result.add_options()("updateSdpThresholdMax",
                              boost::program_options::value<El::BigFloat>(&update_sdp_threshold_max)->default_value(1),
                              "Take a step in the external parameters, "
                              "that is to regenerate the sdp files if the step size is smaller than the threshold. "
                              "The default value is set to 1.");  
  result.add_options()("updateSdpThresholdMin",
                              boost::program_options::value<El::BigFloat>(&update_sdp_threshold_min)->default_value(0),
                              "Take a step in the external parameters, "
                              "that is to regenerate the sdp files if the step size is bigger than the threshold. "
                              "The default value is set to 0.");  
  result.add_options()("updateAtDualityGap",
                              boost::program_options::value<El::BigFloat>(&updateSDP_dualityGapThreshold)->default_value(0),
                              "If updateAtDualityGap is setted to >0, the solver will run until duality<updateAtDualityGap.");
  result.add_options()("centeringRThreshold",
                              boost::program_options::value<El::BigFloat>(&centeringRThreshold)->default_value(-1),
                              "If positive, run centering steps until R<centeringRThreshold.");

  result.add_options()("dualityGapUpperLimit",
	  boost::program_options::value<El::BigFloat>(&dualityGap_upper_limit)->default_value(0),
	  "If dualityGap_upper_limit is setted to >0, mu will not go beyond dualityGap_upper_limit during climbing."
	  "This should be setted to a proper value if there is a bifurcation in global central path when dualityGap > dualityGap_upper_limit.");

  result.add_options()("findBoundary", 
                              boost::program_options::value<bool>(&find_boundary) -> default_value(0), 
                              "True if the program is required by the user to find a point on the boundary of the island. " 
                              "The default value is set to False.");
  result.add_options()("findBoundaryObjThreshold",
                              boost::program_options::value<El::BigFloat>(&find_boundary_obj_threshold)->default_value(0),
                              "Continue to move towards the boundary if the primal and dual objectives are not sufficiently close to zero. " 
                              "The default value is set to 0.");
  result.add_options()("searchDirection",
                              boost::program_options::value<std::vector<El::BigFloat>>(&search_direction)->multitoken(), 
                              "User-specified directional vector in which the program will looks for a zero. ");


  result.add_options()("useExactHessian",
                              boost::program_options::value<bool>(&use_exact_hessian),
                              "To reinitialize the BFGS hessian with the exact one. ");
  result.add_options()("prevGradientBFGS",
                              boost::program_options::value<std::vector<El::BigFloat>>(&prev_grad)->multitoken(),
                              "The gradient of the Largrangian in the last iteration, used to update the hessian. ");  
  result.add_options()("prevExternalStep",
                              boost::program_options::value<std::vector<El::BigFloat>>(&prev_step)->multitoken(),
                              "the step taken by the last iteration, used to update the hessian. "); 
  result.add_options()("prevHessianBFGS",
                              boost::program_options::value<std::vector<El::BigFloat>>(&hess_BFGS)->multitoken(),
                              "Hessian approximated by BFGS. "); 
  result.add_options()("useHmixedforBFGS",
                              boost::program_options::value<bool>(&use_Hmixed_for_BFGS)->default_value(false),
                              "if True, hess_mixed will be used for hess_BFGS");



  result.add_options()("boundingBoxMax",
                              boost::program_options::value<std::vector<El::BigFloat>>(&bounding_box_max)->multitoken(),
                              "The upper bound of the external variables.");
  result.add_options()("boundingBoxMin",
                              boost::program_options::value<std::vector<El::BigFloat>>(&bounding_box_min)->multitoken(),
                              "The lower bound of the external variables.");
  result.add_options()("fixExtParamDirection",
                    	      boost::program_options::value<bool>(&fix_ext_param_direction)->default_value(0),
	                      "True if the program dp is given by --searchDirection. ");

  result.add_options()("printMore",
	                      boost::program_options::bool_switch(&printMore) ->default_value(true),
	                      "Print R error for each iteration");

  result.add_options()("oldSchurDir",
	  boost::program_options::value<boost::filesystem::path>(&old_schur_path),
	  "Directory containing Schur complement of the previous Newtonian step.");

  result.add_options()("oldSDPDir",
	  boost::program_options::value<boost::filesystem::path>(&old_sdp_path),
	  "Directory containing SDP of the previous Newtonian step.");


  // new strategy parameters
  result.add_options()("betaScanMin",
	  boost::program_options::value<El::BigFloat>(&beta_scan_min)->default_value(0.1),
	  "beta scan min value.");
  result.add_options()("betaScanMax",
	  boost::program_options::value<El::BigFloat>(&beta_scan_max)->default_value(1.01),
	  "beta scan max value.");
  result.add_options()("betaScanStep",
	  boost::program_options::value<El::BigFloat>(&beta_scan_step)->default_value(0.1),
	  "beta scan step.");
  result.add_options()("stepMinThreshold",
	  boost::program_options::value<El::BigFloat>(&step_min_threshold)->default_value(0.1),
	  "step size min.");
  result.add_options()("stepMaxThreshold",
	  boost::program_options::value<El::BigFloat>(&step_max_threshold)->default_value(0.6),
	  "step size max.");
  result.add_options()("muIShift",
	  boost::program_options::value<El::BigFloat>(&lagrangian_muI_shift)->default_value(0.2),
	  "mu*I shift in Lagrangian.");
  result.add_options()("maxClimbing",
	  boost::program_options::value<int>(&max_climbing)->default_value(1),
	  "max climbing allowed.");
  result.add_options()("betaClimbing",
	  boost::program_options::value<El::BigFloat>(&beta_climbing)->default_value(2),
	  "beta value for climbing.");
  result.add_options()("betamuLogDetX",
	  boost::program_options::value<El::BigFloat>(&beta_for_mu_logdetX)->default_value(-1),
	  "beta for mu*log(det(X)) term. if beta=-1, the solver takes value from beta scan.");
  result.add_options()("gradientWithLogDetX",
	  boost::program_options::value<bool>(&gradientWithLogDetX)->default_value(false),
	  "if gradientWithLogDetX== true, the gradient N_p is computed with the mu*log(det(X)) term.");

  result.add_options()("finiteDualityGapTarget",
	  boost::program_options::value<El::BigFloat>(&finite_dGap_target)->default_value(-1),
	  "if finiteDualityGapTarget>0, the solver will ignore beta scan and try to move to the specified finite duality gap.");

  result.add_options()("BFGSPartialUpdate",
	  boost::program_options::value<El::BigFloat>(&BFGS_partial_update_reduction)->default_value(-1),
	  "if BFGSPartialUpdate>0, the solver will update BFGS Hessian paritially when the full update is non-positive.");

  result.add_options()("navigatorValueShift",
	  boost::program_options::value<El::BigFloat>(&navigatorValueShift)->default_value(0),
	  "The navigator function will be shifted by this value. This only affects the findBoundary=True run.");


  result.add(solver_parameters.options()); 

  return result;
}
