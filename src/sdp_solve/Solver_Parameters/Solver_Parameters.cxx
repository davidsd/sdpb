#include "../Solver_Parameters.hxx"

#include "Memory_Limit.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <boost/program_options.hpp>

namespace fs = std::filesystem;

Solver_Parameters::Solver_Parameters(const Environment &env)
    : memory_limit_translator(env.initial_node_mem_available())
{}

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
    "preprocessing the input PMP files with 'pmp2sdp'.  GMP will round "
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
    "maxMemory",
    boost::program_options::value<std::string>()
      ->notifier([this](const std::string &s) {
        this->max_memory = memory_limit_translator.from_string(s);
      })
      ->default_value("80%"),
    "Maximum amount of memory that can be used for SDPB allocations, "
    "in bytes."
    " Optional suffixes: B (bytes), K/KB/KiB (kilobytes), M/MB/MiB "
    "(megabytes), "
    "G/GB/GiB (gigabytes), % (percents of MemAvailable at program start)");
  result.add_options()(
    "maxSharedMemory",
    boost::program_options::value<std::string>()
      ->notifier([this](const std::string &s) {
        this->max_shared_memory = memory_limit_translator.from_string(s);
      })
      ->default_value("0"),
    "Maximum amount of memory that can be used for MPI shared windows, "
    "in bytes."
    " Optional suffixes: B (bytes), K/KB/KiB (kilobytes), M/MB/MiB "
    "(megabytes), "
    "G/GB/GiB (gigabytes), % (percents of MemAvailable at program start)");
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
    boost::program_options::value<fs::path>(&checkpoint_out),
    "Checkpoints are saved to this directory every checkpointInterval.  Set "
    "to the empty string to inhibit checkpoints.  Defaults to sdpDir with "
    "'.ck' extension.");
  result.add_options()("initialCheckpointDir,i",
                       boost::program_options::value<fs::path>(&checkpoint_in),
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

boost::property_tree::ptree to_property_tree(const Solver_Parameters &p)
{
  boost::property_tree::ptree result;

  result.put("maxIterations", p.max_iterations);
  result.put("maxRuntime", p.max_runtime);
  result.put("maxMemory", p.max_memory, p.memory_limit_translator);
  result.put("maxSharedMemory", p.max_shared_memory,
             p.memory_limit_translator);
  result.put("checkpointInterval", p.checkpoint_interval);
  result.put("findPrimalFeasible", p.find_primal_feasible);
  result.put("findDualFeasible", p.find_dual_feasible);
  result.put("detectPrimalFeasibleJump", p.detect_primal_feasible_jump);
  result.put("detectDualFeasibleJump", p.detect_dual_feasible_jump);
  result.put("precision", p.precision);
  result.put("precision_actual", mpf_get_default_prec());
  result.put("dualityGapThreshold", p.duality_gap_threshold);
  result.put("primalErrorThreshold", p.primal_error_threshold);
  result.put("dualErrorThreshold", p.dual_error_threshold);
  result.put("initialMatrixScalePrimal", p.initial_matrix_scale_primal);
  result.put("initialMatrixScaleDual", p.initial_matrix_scale_dual);
  result.put("feasibleCenteringParameter", p.feasible_centering_parameter);
  result.put("infeasibleCenteringParameter", p.infeasible_centering_parameter);
  result.put("stepLengthReduction", p.step_length_reduction);
  result.put("maxComplementarity", p.max_complementarity);
  result.put("initialCheckpointDir", p.checkpoint_in.string());
  result.put("checkpointDir", p.checkpoint_out.string());

  return result;
}

std::ostream &operator<<(std::ostream &os, const Solver_Parameters &p)
{
  os << std::boolalpha << "maxIterations                = " << p.max_iterations
     << '\n'
     << "maxRuntime                   = " << p.max_runtime << '\n'
     << "checkpointInterval           = " << p.checkpoint_interval << '\n'
     << "maxMemory                    = " << p.max_memory << '\n'
     << "maxSharedMemory              = " << p.max_shared_memory << '\n'
     << "findPrimalFeasible           = " << p.find_primal_feasible << '\n'
     << "findDualFeasible             = " << p.find_dual_feasible << '\n'
     << "detectPrimalFeasibleJump     = " << p.detect_primal_feasible_jump
     << '\n'
     << "detectDualFeasibleJump       = " << p.detect_dual_feasible_jump
     << '\n'
     << "precision(actual)            = " << p.precision << "("
     << mpf_get_default_prec() << ")" << '\n'

     << "dualityGapThreshold          = " << p.duality_gap_threshold << '\n'
     << "primalErrorThreshold         = " << p.primal_error_threshold << '\n'
     << "dualErrorThreshold           = " << p.dual_error_threshold << '\n'
     << "initialMatrixScalePrimal     = " << p.initial_matrix_scale_primal
     << '\n'
     << "initialMatrixScaleDual       = " << p.initial_matrix_scale_dual
     << '\n'
     << "feasibleCenteringParameter   = " << p.feasible_centering_parameter
     << '\n'
     << "infeasibleCenteringParameter = " << p.infeasible_centering_parameter
     << '\n'
     << "stepLengthReduction          = " << p.step_length_reduction << '\n'
     << "maxComplementarity           = " << p.max_complementarity << '\n'
     << "initialCheckpointDir         = " << p.checkpoint_in << '\n'
     << "checkpointDir                = " << p.checkpoint_out << '\n';
  return os;
}
