#include "../Dynamical_Solver_Parameters.hxx"

#include <boost/program_options.hpp>
#include <vector>

boost::program_options::options_description
Dynamical_Solver_Parameters::options()
{
  boost::program_options::options_description result(
    "Dynamic Solver parameters");
  // We set default parameters using El::BigFloat("1e-10",10)
  // rather than a straight double precision 1e-10 so that results
  // are more reproducible at high precision.  Using double
  // precision defaults results in differences of about 1e-15 in
  // primalObjective after one step.

  result.add_options()(
    "newSdpDirs",
    boost::program_options::value<std::filesystem::path>(&new_sdp_path)
      ->required(),
    "Directory containing the preprocessed SDP data files around the center "
    "SDP in external parameter space.");

  result.add_options()(
    "externalParamInfinitestimal",
    boost::program_options::value<El::BigFloat>(&alpha)->default_value(1),
    "Step size in the external-parameter space to generate the new SDP data "
    "files. "
    "The default value is set to 1.");
  result.add_options()(
    "numExternalParams",
    boost::program_options::value<int>(&n_external_parameters)
      ->default_value(0),
    "The number of external parameters to be varied in each iteration of the "
    "dynamical SDP. "
    "The default value is set to 0.");
  result.add_options()(
    "totalIterationCount",
    boost::program_options::value<size_t>(&total_iterations)->default_value(0),
    "The number of total iterations finished before this solver starts to "
    "run. "
    "The default value is set to 0.");
  result.add_options()(
    "updateSdpThresholdMax",
    boost::program_options::value<El::BigFloat>(&update_sdp_threshold_max)
      ->default_value(1),
    "Take a step in the external parameters, "
    "that is to regenerate the sdp files if the step size is smaller than the "
    "threshold. "
    "The default value is set to 1.");
  result.add_options()(
    "updateSdpThresholdMin",
    boost::program_options::value<El::BigFloat>(&update_sdp_threshold_min)
      ->default_value(0),
    "Take a step in the external parameters, "
    "that is to regenerate the sdp files if the step size is bigger than the "
    "threshold. "
    "The default value is set to 0.");
  result.add_options()(
    "centeringRThreshold",
    boost::program_options::value<El::BigFloat>(&centeringRThreshold)
      ->default_value(-1),
    "If positive, run centering steps until R<centeringRThreshold.");

  result.add_options()(
    "dualityGapUpperLimit",
    boost::program_options::value<El::BigFloat>(&dualityGap_upper_limit)
      ->default_value(0),
    "If dualityGap_upper_limit is setted to >0, mu will not go beyond "
    "dualityGap_upper_limit during climbing."
    "This should be setted to a proper value if there is a bifurcation in "
    "global central path when dualityGap > dualityGap_upper_limit.");

  result.add_options()(
    "findBoundaryObjThreshold",
    boost::program_options::value<El::BigFloat>(&find_boundary_obj_threshold)
      ->default_value(0),
    "Continue to move towards the boundary if the primal and dual objectives "
    "are not sufficiently close to zero. "
    "The default value is set to 0.");
  result.add_options()(
    "findBoundaryDirection",
    boost::program_options::value<std::vector<El::BigFloat>>(&search_direction)
      ->multitoken(),
    "User-specified directional vector in which the program will look for a "
    "zero. "
    "Without specify this option, the program will look for the minimum of "
    "the navigator function.");

  result.add_options()(
    "useExactHessian", boost::program_options::value<bool>(&use_exact_hessian),
    "To reinitialize the BFGS hessian with the exact one. ");
  result.add_options()(
    "prevGradientBFGS",
    boost::program_options::value<std::vector<El::BigFloat>>(&prev_grad)
      ->multitoken(),
    "The gradient of the Largrangian in the last iteration, used to update "
    "the hessian. ");
  result.add_options()(
    "prevExternalStep",
    boost::program_options::value<std::vector<El::BigFloat>>(&prev_step)
      ->multitoken(),
    "the step taken by the last iteration, used to update the hessian. ");
  result.add_options()(
    "prevHessianBFGS",
    boost::program_options::value<std::vector<El::BigFloat>>(&hess_BFGS)
      ->multitoken(),
    "Hessian approximated by BFGS. ");
  result.add_options()(
    "prevHessianBFGSpp",
    boost::program_options::value<std::vector<El::BigFloat>>(&hess_BFGS_pp)
      ->multitoken(),
    "H_pp approximated by BFGS. ");

  result.add_options()(
    "boundingBoxMax",
    boost::program_options::value<std::vector<El::BigFloat>>(&bounding_box_max)
      ->multitoken(),
    "The upper bound of the external variables.");
  result.add_options()(
    "boundingBoxMin",
    boost::program_options::value<std::vector<El::BigFloat>>(&bounding_box_min)
      ->multitoken(),
    "The lower bound of the external variables.");
  result.add_options()(
    "fixExtParamDirection",
    boost::program_options::value<bool>(&fix_ext_param_direction)
      ->default_value(0),
    "True if the program dp is given by --searchDirection. ");

  result.add_options()(
    "printMore",
    boost::program_options::bool_switch(&printMore)->default_value(true),
    "Print R error for each iteration");

  result.add_options()(
    "oldSchurDir",
    boost::program_options::value<std::filesystem::path>(&old_schur_path),
    "Directory containing Schur complement of the previous Newtonian step.");

  result.add_options()(
    "oldSDPDir",
    boost::program_options::value<std::filesystem::path>(&old_sdp_path),
    "Directory containing SDP of the previous Newtonian step.");

  // new strategy parameters
  result.add_options()(
    "betaScanMin",
    boost::program_options::value<El::BigFloat>(&beta_scan_min)
      ->default_value(0.1),
    "beta scan min value.");
  result.add_options()(
    "betaScanMax",
    boost::program_options::value<El::BigFloat>(&beta_scan_max)
      ->default_value(1.01),
    "beta scan max value.");
  result.add_options()(
    "betaScanStep",
    boost::program_options::value<El::BigFloat>(&beta_scan_step)
      ->default_value(0.1),
    "beta scan step.");
  result.add_options()(
    "stepMinThreshold",
    boost::program_options::value<El::BigFloat>(&step_min_threshold)
      ->default_value(0.1),
    "step size min.");
  result.add_options()(
    "stepMaxThreshold",
    boost::program_options::value<El::BigFloat>(&step_max_threshold)
      ->default_value(0.6),
    "step size max.");
  result.add_options()(
    "primalDualObjWeight",
    boost::program_options::value<El::BigFloat>(&lagrangian_muI_shift)
      ->default_value(0.2),
    "mu*I shift in Lagrangian.");
  result.add_options()(
    "maxClimbingSteps",
    boost::program_options::value<int>(&max_climbing)->default_value(1),
    "max climbing allowed.");
  result.add_options()(
    "betaClimbing",
    boost::program_options::value<El::BigFloat>(&beta_climbing)
      ->default_value(2),
    "beta value for climbing.");

  result.add_options()(
    "navigatorWithLogDetX",
    boost::program_options::value<bool>(&navigatorWithLogDetX)
      ->default_value(false),
    "if gradientWithLogDetX== true, the navigator value is computed with the "
    "mu*log(det(X)) term.");
  result.add_options()(
    "gradientWithLogDetX",
    boost::program_options::value<bool>(&gradientWithLogDetX)
      ->default_value(true),
    "if gradientWithLogDetX== true, the gradient N_p is computed with the "
    "mu*log(det(X)) term.");

  result.add_options()(
    "finiteDualityGapTarget",
    boost::program_options::value<El::BigFloat>(&finite_dGap_target)
      ->default_value(-1),
    "if finiteDualityGapTarget>0, the solver will ignore beta scan and try to "
    "move to the specified finite duality gap.");

  result.add_options()(
    "BFGSPartialUpdate",
    boost::program_options::value<El::BigFloat>(&BFGS_partial_update_reduction)
      ->default_value(-1),
    "if BFGSPartialUpdate>0, the solver will update BFGS Hessian paritially "
    "when the full update is non-positive.");

  result.add_options()(
    "navigatorValueShift",
    boost::program_options::value<El::BigFloat>(&navigatorValueShift)
      ->default_value(0),
    "Experimental : The navigator function will be shifted by this value. "
    "This only affects the findBoundary=True run.");

  result.add_options()(
    "navigatorAutomaticShift",
    boost::program_options::value<bool>(&navigatorAutomaticShiftQ)
      ->default_value(false),
    "Experimental : shift the navigator by dGap*log(dGap/dim(X)).");

  result.add_options()(
    "stickToGCP",
    boost::program_options::value<bool>(&stickToGCPQ)->default_value(false),
    "Experimental : stick to GCP.");

  result.add_options()(
    "optimalbeta",
    boost::program_options::value<bool>(&optimalbetaQ)->default_value(false),
    "Experimental : optimal beta.");

  result.add_options()(
    "climbingRecomputeExtParam",
    boost::program_options::value<bool>(&climbingRecomputeExtParamQ)
      ->default_value(true),
    "Experimental : recompute external parameter after clibming. Default : "
    "true .");

  result.add_options()(
    "returnCheckpointOnLCP",
    boost::program_options::value<bool>(&returnCheckpointOnLCP)
      ->default_value(false),
    "Experimental : the code will return checkpoint on local central path.");

  result.add(solver_parameters.options());

  return result;
}

namespace po = boost::program_options;

PrecParameters::PrecParameters(int argc, char *argv[])
{
  po::options_description cmd_line_options("precision option");

  cmd_line_options.add_options()(
    "precision",
    boost::program_options::value<size_t>(&prec)->default_value(768),
    "The precision, in the number of bits, for numbers in the "
    "computation. "
    " This should be less than or equal to the precision used when "
    "preprocessing the XML input files with 'pvm2sdp'.  GMP will round "
    "this up to a multiple of 32 or 64, depending on the system.");

  po::variables_map variables_map;
  try
    {
      //po::store(po::parse_command_line(argc, argv, cmd_line_options, po::command_line_style::unix_style ^ po::command_line_style::allow_short),
      //	variables_map);

      po::store(po::command_line_parser(argc, argv)
                  .options(cmd_line_options)
                  .allow_unregistered()
                  .run(),
                variables_map);

      po::notify(variables_map);
    }
  catch(po::error &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
