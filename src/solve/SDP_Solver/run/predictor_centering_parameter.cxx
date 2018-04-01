#include "../../SDP_Solver.hxx"

// Centering parameter \beta_p for the predictor step
Real predictor_centering_parameter(const SDP_Solver_Parameters &parameters,
                                   const bool is_primal_dual_feasible)
{
  return is_primal_dual_feasible ? Real(0)
                                 : parameters.infeasible_centering_parameter;
}
