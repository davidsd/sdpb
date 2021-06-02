#include "../Solver_Parameters.hxx"

#include <boost/program_options.hpp>

boost::program_options::options_description Solver_Parameters::options()
{
  boost::program_options::options_description result("Solver parameters");
  // We set default parameters using El::BigFloat("1e-10",10)
  // rather than a straight double precision 1e-10 so that results
  // are more reproducible at high precision.  Using double
  // precision defaults results in differences of about 1e-15 in
  // primalObjective after one step.

  result.add_options()(
    "precision",
    boost::program_options::value<size_t>(&precision)->default_value(400),
    "The precision, in the number of bits, for numbers in the "
    "computation. "
    " This should be less than or equal to the precision used when "
    "preprocessing the XML input files with 'pvm2sdp'.  GMP will round "
    "this up to a multiple of 32 or 64, depending on the system.");
  result.add_options()(
    "findPrimalFeasible",
    boost::program_options::bool_switch(&find_primal_feasible)
      ->default_value(false),
    "Terminate once a primal feasible solution is found.");
  result.add_options()("findDualFeasible",
                       boost::program_options::bool_switch(&find_dual_feasible)
                         ->default_value(false),
                       "Terminate once a dual feasible solution is found.");
  result.add_options()(
    "detectPrimalFeasibleJump",
    boost::program_options::bool_switch(&detect_primal_feasible_jump)
      ->default_value(false),
    "Terminate if a primal-step of 1 is taken. This often indicates that "
    "a "
    "primal feasible solution would be found if the precision were high "
    "enough. Try increasing either primalErrorThreshold or precision "
    "and run from the latest checkpoint.");
  result.add_options()(
    "detectDualFeasibleJump",
    boost::program_options::bool_switch(&detect_dual_feasible_jump)
      ->default_value(false),
    "Terminate if a dual-step of 1 is taken. This often indicates that a "
    "dual feasible solution would be found if the precision were high "
    "enough. Try increasing either dualErrorThreshold or precision "
    "and run from the latest checkpoint.");
  result.add_options()(
    "maxIterations",
    boost::program_options::value<int64_t>(&max_iterations)->default_value(500),
    "Maximum number of iterations to run the solver.");
  result.add_options()("maxRuntime",
                       boost::program_options::value<int64_t>(&max_runtime)
                         ->default_value(std::numeric_limits<int64_t>::max()),
                       "Maximum amount of time to run the solver in seconds.");
  result.add_options()(
    "dualityGapThreshold",
    boost::program_options::value<El::BigFloat>(&duality_gap_threshold)
      ->default_value(El::BigFloat("1e-30", 10)),
    "Threshold for duality gap (roughly the difference in primal and dual "
    "objective) at which the solution is considered "
    "optimal. Corresponds to SDPA's epsilonStar.");
  result.add_options()(
    "primalErrorThreshold",
    boost::program_options::value<El::BigFloat>(&primal_error_threshold)
      ->default_value(El::BigFloat("1e-30", 10)),
    "Threshold for feasibility of the primal problem. Corresponds to "
    "SDPA's epsilonBar.");
  result.add_options()(
    "dualErrorThreshold",
    boost::program_options::value<El::BigFloat>(&dual_error_threshold)
      ->default_value(El::BigFloat("1e-30", 10)),
    "Threshold for feasibility of the dual problem. Corresponds to SDPA's "
    "epsilonBar.");
  result.add_options()(
    "initialMatrixScalePrimal",
    boost::program_options::value<El::BigFloat>(&initial_matrix_scale_primal)
      ->default_value(El::BigFloat("1e20", 10)),
    "The primal matrix X begins at initialMatrixScalePrimal times the "
    "identity matrix. Corresponds to SDPA's lambdaStar.");
  result.add_options()(
    "initialMatrixScaleDual",
    boost::program_options::value<El::BigFloat>(&initial_matrix_scale_dual)
      ->default_value(El::BigFloat("1e20", 10)),
    "The dual matrix Y begins at initialMatrixScaleDual times the "
    "identity matrix. Corresponds to SDPA's lambdaStar.");
  result.add_options()(
    "feasibleCenteringParameter",
    boost::program_options::value<El::BigFloat>(&feasible_centering_parameter)
      ->default_value(El::BigFloat("0.1", 10)),
    "Shrink the complementarity X Y by this factor when the primal and "
    "dual "
    "problems are feasible. Corresponds to SDPA's betaStar.");
  result.add_options()(
    "infeasibleCenteringParameter",
    boost::program_options::value<El::BigFloat>(
      &infeasible_centering_parameter)
      ->default_value(El::BigFloat("0.3", 10)),
    "Shrink the complementarity X Y by this factor when either the primal "
    "or dual problems are infeasible. Corresponds to SDPA's betaBar.");
  result.add_options()(
    "stepLengthReduction",
    boost::program_options::value<El::BigFloat>(&step_length_reduction)
      ->default_value(El::BigFloat("0.7", 10)),
    "Shrink each newton step by this factor (smaller means slower, more "
    "stable convergence). Corresponds to SDPA's gammaStar.");
  result.add_options()(
    "minPrimalStep",
    boost::program_options::value<El::BigFloat>(&min_primal_step)
      ->default_value(El::BigFloat(0)),
    "Terminate if the primal step size becomes smaller than this value.");
  result.add_options()(
    "minDualStep",
    boost::program_options::value<El::BigFloat>(&min_dual_step)
      ->default_value(El::BigFloat(0)),
    "Terminate if the dual step size becomes smaller than this value.");
  result.add_options()(
    "maxComplementarity",
    boost::program_options::value<El::BigFloat>(&max_complementarity)
      ->default_value(El::BigFloat("1e100", 10)),
    "Terminate if the complementarity mu = Tr(X Y)/dim(X) "
    "exceeds this value.");
  result.add_options()(
    "checkpointDir,c",
    boost::program_options::value<boost::filesystem::path>(&checkpoint_out),
    "Checkpoints are saved to this directory every checkpointInterval.  Set "
    "to the empty string to inhibit checkpoints.  Defaults to sdpDir with "
    "'.ck' extension.");
  result.add_options()(
    "initialCheckpointDir,i",
    boost::program_options::value<boost::filesystem::path>(&checkpoint_in),
    "The initial checkpoint directory to load. Defaults to "
    "checkpointDir.");
  result.add_options()(
    "checkpointInterval",
    boost::program_options::value<int64_t>(&checkpoint_interval)
      ->default_value(3600),
    "Save checkpoints to checkpointDir every checkpointInterval "
    "seconds.");

  return result;
}
