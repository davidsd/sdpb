#include "../../SDP_Solver_Terminate_Reason.hxx"
#include "../../../SDP_Solver_Parameters.hxx"

void compute_feasible_and_termination(
  const SDP_Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration, bool &is_primal_and_dual_feasible,
  SDP_Solver_Terminate_Reason &result, bool &terminate_now)
{
  const bool is_primal_feasible(primal_error
                                < parameters.primal_error_threshold),
    is_dual_feasible(dual_error < parameters.dual_error_threshold),
    is_optimal(duality_gap < parameters.duality_gap_threshold);

  terminate_now=true;
  if(is_primal_feasible && is_dual_feasible && is_optimal)
    {
      result = SDP_Solver_Terminate_Reason::PrimalDualOptimal;
    }
  else if(is_primal_feasible && parameters.find_primal_feasible)
    {
      result = SDP_Solver_Terminate_Reason::PrimalFeasible;
    }
  else if(is_dual_feasible && parameters.find_dual_feasible)
    {
      result = SDP_Solver_Terminate_Reason::DualFeasible;
    }
  else if(primal_step_length == El::BigFloat(1)
          && parameters.detect_primal_feasible_jump)
    {
      result = SDP_Solver_Terminate_Reason::PrimalFeasibleJumpDetected;
    }
  else if(dual_step_length == El::BigFloat(1)
          && parameters.detect_dual_feasible_jump)
    {
      result = SDP_Solver_Terminate_Reason::DualFeasibleJumpDetected;
    }
  else if(iteration > parameters.max_iterations)
    {
      result = SDP_Solver_Terminate_Reason::MaxIterationsExceeded;
    }
  else
    {
      terminate_now=false;
    }

  is_primal_and_dual_feasible = (is_primal_feasible && is_dual_feasible);
}
