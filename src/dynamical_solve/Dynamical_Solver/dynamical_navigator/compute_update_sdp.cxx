#include "../../../sdp_solve/Solver_Parameters.hxx"
#include "../../Dynamical_Solver_Parameters.hxx"

void compute_update_sdp(
  const Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  bool &update_sdp)
{
  update_sdp = true; 

}

