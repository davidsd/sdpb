#pragma once
// Parameters to control the behavior of an SDPSolver. See the manual
// for a detailed description of each.
//

#include "types.hxx"

#include <iostream>

class SDP_Solver_Parameters
{
public:
  int maxIterations;
  int maxRuntime;
  int checkpointInterval;
  bool noFinalCheckpoint;
  bool findPrimalFeasible;
  bool findDualFeasible;
  bool detectPrimalFeasibleJump;
  bool detectDualFeasibleJump;
  int precision;
  int maxThreads;
  Real dualityGapThreshold;
  Real primalErrorThreshold;
  Real dualErrorThreshold;
  Real initialMatrixScalePrimal;
  Real initialMatrixScaleDual;
  Real feasibleCenteringParameter;
  Real infeasibleCenteringParameter;
  Real stepLengthReduction;
  Real choleskyStabilizeThreshold;
  Real maxComplementarity;
  // Set the precision of all Real parameters to equal 'precision'.
  // This is necessary because 'precision' might be set (via the
  // command line or a file) after initializing other parameters.
  //
  void resetPrecision()
  {
    dualityGapThreshold.set_prec(precision);
    primalErrorThreshold.set_prec(precision);
    dualErrorThreshold.set_prec(precision);
    initialMatrixScalePrimal.set_prec(precision);
    initialMatrixScaleDual.set_prec(precision);
    feasibleCenteringParameter.set_prec(precision);
    infeasibleCenteringParameter.set_prec(precision);
    stepLengthReduction.set_prec(precision);
    choleskyStabilizeThreshold.set_prec(precision);
    maxComplementarity.set_prec(precision);
  }

  friend std::ostream &
  operator<<(std::ostream &os, const SDP_Solver_Parameters &p);
};
