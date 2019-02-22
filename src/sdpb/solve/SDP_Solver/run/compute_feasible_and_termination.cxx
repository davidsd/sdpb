#include "../../SDP_Solver_Terminate_Reason.hxx"
#include "../../../SDP_Solver_Parameters.hxx"

void compute_feasible_and_termination(
  const SDP_Solver_Parameters &parameters,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_dual_feasible, bool &is_optimal,
  SDP_Solver_Terminate_Reason &result, bool &terminate_now)
{
  is_dual_feasible = dual_error < parameters.dual_error_threshold;
  is_optimal = duality_gap < parameters.duality_gap_threshold;

  terminate_now = true;
  if(is_dual_feasible && parameters.find_dual_feasible)
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
  else if(std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::high_resolution_clock::now() - solver_start_time)
            .count()
          >= parameters.max_runtime)
    {
      result = SDP_Solver_Terminate_Reason::MaxRuntimeExceeded;
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
