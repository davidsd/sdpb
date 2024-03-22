#include "../Dynamical_Solver_Terminate_Reason.hxx"
#include "sdpb_util/assert.hxx"

std::ostream &
operator<<(std::ostream &os, const Dynamical_Solver_Terminate_Reason &r)
{
  switch(r)
    {
    case Dynamical_Solver_Terminate_Reason::PrimalDualOptimal:
      os << "found primal-dual optimal solution";
      break;
    case Dynamical_Solver_Terminate_Reason::PrimalFeasible:
      os << "found primal feasible solution";
      break;
    case Dynamical_Solver_Terminate_Reason::DualFeasible:
      os << "found dual feasible solution";
      break;
    case Dynamical_Solver_Terminate_Reason::PrimalFeasibleJumpDetected:
      os << "primal feasible jump detected";
      break;
    case Dynamical_Solver_Terminate_Reason::DualFeasibleJumpDetected:
      os << "dual feasible jump detected";
      break;
    case Dynamical_Solver_Terminate_Reason::MaxIterationsExceeded:
      os << "maxIterations exceeded";
      break;
    case Dynamical_Solver_Terminate_Reason::MaxRuntimeExceeded:
      os << "maxRuntime exceeded";
      break;
    case Dynamical_Solver_Terminate_Reason::MaxComplementarityExceeded:
      os << "maxComplementarity exceeded";
      break;
    case Dynamical_Solver_Terminate_Reason::PrimalStepTooSmall:
      os << "primal step too small";
      break;
    case Dynamical_Solver_Terminate_Reason::DualStepTooSmall:
      os << "dual step too small";
      break;
    case Dynamical_Solver_Terminate_Reason::UpdateSDPs:
      os << "met criteria to update sdp input files";
      break;
    case Dynamical_Solver_Terminate_Reason::SIGTERM_Received:
      os << "SIGTERM signal received";
      break;
    default: RUNTIME_ERROR("Unknown SDP_Solver_Terminate_Reason=", r);
    }
  return os;
}
