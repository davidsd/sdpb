#include "../../SDP_Solver_Terminate_Reason.hxx"
#include "../../../SDP_Solver_Parameters.hxx"

void compute_feasible_and_termination(
  const SDP_Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  SDP_Solver_Terminate_Reason &terminate_reason, bool &terminate_now)
{
  const bool is_dual_feasible(dual_error < parameters.dual_error_threshold),
    is_primal_feasible(primal_error < parameters.primal_error_threshold);
  is_primal_and_dual_feasible = (is_primal_feasible && is_dual_feasible);

  const bool is_optimal(duality_gap < parameters.duality_gap_threshold);

  terminate_now = true;
  if(is_primal_and_dual_feasible && is_optimal)
    {
      terminate_reason = SDP_Solver_Terminate_Reason::PrimalDualOptimal;
    }
  else if(is_dual_feasible && parameters.find_dual_feasible)
    {
      terminate_reason = SDP_Solver_Terminate_Reason::DualFeasible;
    }
  else if(is_primal_feasible && parameters.find_primal_feasible)
    {
      terminate_reason = SDP_Solver_Terminate_Reason::PrimalFeasible;
    }
  else if(dual_step_length == El::BigFloat(1)
          && parameters.detect_dual_feasible_jump)
    {
      terminate_reason = SDP_Solver_Terminate_Reason::DualFeasibleJumpDetected;
    }
  else if(primal_step_length == El::BigFloat(1)
          && parameters.detect_primal_feasible_jump)
    {
      terminate_reason
        = SDP_Solver_Terminate_Reason::PrimalFeasibleJumpDetected;
    }
  else if(iteration > parameters.max_iterations)
    {
      terminate_reason = SDP_Solver_Terminate_Reason::MaxIterationsExceeded;
    }
  else if(std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::high_resolution_clock::now() - solver_start_time)
            .count()
          >= parameters.max_runtime)
    {
      terminate_reason = SDP_Solver_Terminate_Reason::MaxRuntimeExceeded;
    }
  else
    {
      terminate_now = false;
    }
  El::byte terminate_byte(terminate_now);
  // Time varies between cores, so follow the decision of the root.
  El::mpi::Broadcast(terminate_byte, 0, El::mpi::COMM_WORLD);
  terminate_now = static_cast<bool>(terminate_byte);
}
