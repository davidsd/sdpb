#include "../SDP_Solver_Terminate_Reason.hxx"

std::ostream &
operator<<(std::ostream &os, const SDP_Solver_Terminate_Reason &r)
{
  switch(r)
    {
    case SDP_Solver_Terminate_Reason::PrimalDualOptimal:
      os << "found primal-dual optimal solution";
      break;
    case SDP_Solver_Terminate_Reason::PrimalFeasible:
      os << "found primal feasible solution";
      break;
    case SDP_Solver_Terminate_Reason::DualFeasible:
      os << "found dual feasible solution";
      break;
    case SDP_Solver_Terminate_Reason::PrimalFeasibleJumpDetected:
      os << "primal feasible jump detected";
      break;
    case SDP_Solver_Terminate_Reason::DualFeasibleJumpDetected:
      os << "dual feasible jump detected";
      break;
    case SDP_Solver_Terminate_Reason::MaxIterationsExceeded:
      os << "maxIterations exceeded";
      break;
    case SDP_Solver_Terminate_Reason::MaxRuntimeExceeded:
      os << "maxRuntime exceeded";
      break;
    case SDP_Solver_Terminate_Reason::MaxComplementarityExceeded:
      os << "maxComplementarity exceeded";
      break;
    case SDP_Solver_Terminate_Reason::PrimalStepTooSmall:
      os << "primal step too small";
      break;
    case SDP_Solver_Terminate_Reason::DualStepTooSmall:
      os << "dual step too small";
      break;
    }
  return os;
}
