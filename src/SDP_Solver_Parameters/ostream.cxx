#include "../SDP_Solver_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const SDP_Solver_Parameters &p)
{
  os << std::boolalpha;
  os << "maxIterations                = " << p.maxIterations << '\n';
  os << "maxRuntime                   = " << p.maxRuntime << '\n';
  os << "checkpointInterval           = " << p.checkpointInterval << '\n';
  os << "noFinalCheckpoint            = " << p.noFinalCheckpoint << '\n';
  os << "findPrimalFeasible           = " << p.findPrimalFeasible << '\n';
  os << "findDualFeasible             = " << p.findDualFeasible << '\n';
  os << "detectPrimalFeasibleJump     = " << p.detectPrimalFeasibleJump
     << '\n';
  os << "detectDualFeasibleJump       = " << p.detectDualFeasibleJump << '\n';
  os << "precision(actual)            = " << p.precision << "("
     << mpf_get_default_prec() << ")" << '\n';
  os << "dualityGapThreshold          = " << p.dualityGapThreshold << '\n';
  os << "primalErrorThreshold         = " << p.primalErrorThreshold << '\n';
  os << "dualErrorThreshold           = " << p.dualErrorThreshold << '\n';
  os << "initialMatrixScalePrimal     = " << p.initialMatrixScalePrimal
     << '\n';
  os << "initialMatrixScaleDual       = " << p.initialMatrixScaleDual << '\n';
  os << "feasibleCenteringParameter   = " << p.feasibleCenteringParameter
     << '\n';
  os << "infeasibleCenteringParameter = " << p.infeasibleCenteringParameter
     << '\n';
  os << "stepLengthReduction          = " << p.stepLengthReduction << '\n';
  os << "choleskyStabilizeThreshold   = " << p.choleskyStabilizeThreshold
     << '\n';
  os << "maxComplementarity           = " << p.maxComplementarity << '\n';
  return os;
}

