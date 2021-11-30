#pragma once

#include <ostream>

// Reasons for terminating the solver.  See the manual for a detailed
// description of each.

enum class Dynamical_Solver_Terminate_Reason
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
  UpdatingSDPs, 
};

std::ostream &
operator<<(std::ostream &os, const Dynamical_Solver_Terminate_Reason &r);
