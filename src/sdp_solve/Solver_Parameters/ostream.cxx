#include "../Solver_Parameters.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

std::ostream &operator<<(std::ostream &os, const Solver_Parameters &p)
{
  os << std::boolalpha << "maxIterations                = " << p.max_iterations
     << '\n'
     << "maxRuntime                   = " << p.max_runtime << '\n'
     << "checkpointInterval           = " << p.checkpoint_interval << '\n'
     << "maxSharedMemory              = "
     << pretty_print_bytes(p.max_shared_memory_bytes, true) << '\n'
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
