#include "../SDP_Solver_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const SDP_Solver_Parameters &p)
{
  os << std::boolalpha;
  os << "maxIterations                = " << p.max_iterations << '\n';
  os << "maxRuntime                   = " << p.max_runtime << '\n';
  os << "checkpointInterval           = " << p.checkpoint_interval << '\n';
  os << "noFinalCheckpoint            = " << p.no_final_checkpoint << '\n';
  os << "findPrimalFeasible           = " << p.find_primal_feasible << '\n';
  os << "findDualFeasible             = " << p.find_dual_feasible << '\n';
  os << "detectPrimalFeasibleJump     = " << p.detect_primal_feasible_jump
     << '\n';
  os << "detectDualFeasibleJump       = " << p.detect_dual_feasible_jump
     << '\n';
  os << "precision(actual)            = " << p.precision << "("
     << mpf_get_default_prec() << ")" << '\n';
  os << "dualityGapThreshold          = " << p.duality_gap_threshold << '\n';
  os << "primalErrorThreshold         = " << p.primal_error_threshold << '\n';
  os << "dualErrorThreshold           = " << p.dual_error_threshold << '\n';
  os << "initialMatrixScalePrimal     = " << p.initial_matrix_scale_primal
     << '\n';
  os << "initialMatrixScaleDual       = " << p.initial_matrix_scale_dual
     << '\n';
  os << "feasibleCenteringParameter   = " << p.feasible_centering_parameter
     << '\n';
  os << "infeasibleCenteringParameter = " << p.infeasible_centering_parameter
     << '\n';
  os << "stepLengthReduction          = " << p.step_length_reduction << '\n';
  os << "maxComplementarity           = " << p.max_complementarity << '\n';

  os << "dualityGapThreshold          = " << p.duality_gap_threshold_elemental << '\n';
  os << "primalErrorThreshold         = " << p.primal_error_threshold_elemental << '\n';
  os << "dualErrorThreshold           = " << p.dual_error_threshold_elemental << '\n';
  os << "initialMatrixScalePrimal     = " << p.initial_matrix_scale_primal_elemental
     << '\n';
  os << "initialMatrixScaleDual       = " << p.initial_matrix_scale_dual_elemental
     << '\n';
  os << "feasibleCenteringParameter   = " << p.feasible_centering_parameter_elemental
     << '\n';
  os << "infeasibleCenteringParameter = " << p.infeasible_centering_parameter_elemental
     << '\n';
  os << "stepLengthReduction          = " << p.step_length_reduction_elemental << '\n';
  os << "maxComplementarity           = " << p.max_complementarity_elemental << '\n';
  return os;
}
