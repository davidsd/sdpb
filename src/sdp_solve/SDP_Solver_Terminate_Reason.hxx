#pragma once

#include <ostream>

// Reasons for terminating the solver.  See the manual for a detailed
// description of each.

enum class SDP_Solver_Terminate_Reason
{
  PrimalDualOptimal,
  PrimalFeasible,
  DualFeasible,
  PrimalFeasibleJumpDetected,
  DualFeasibleJumpDetected,
  MaxComplementarityExceeded,
  MaxIterationsExceeded,
  MaxRuntimeExceeded,
  PrimalStepTooSmall,
  DualStepTooSmall,
};

std::ostream &
operator<<(std::ostream &os, const SDP_Solver_Terminate_Reason &r);
