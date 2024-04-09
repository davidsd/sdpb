#include "String_To_Bytes_Translator.hxx"
#include "../Solver_Parameters.hxx"

boost::property_tree::ptree to_property_tree(const Solver_Parameters &p)
{
  boost::property_tree::ptree result;

  result.put("maxIterations", p.max_iterations);
  result.put("maxRuntime", p.max_runtime);
  result.put("maxSharedMemory", p.max_shared_memory_bytes,
             String_To_Bytes_Translator());
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
