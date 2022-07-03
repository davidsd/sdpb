#include "../../../sdp_solve/Solver_Parameters.hxx"
#include "../../Dynamical_Solver_Parameters.hxx"

// subroutines to decide whether to update sdps in run_dynamical, before entering dynamical_step
void compute_update_sdp(
  const Dynamical_Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  bool &update_sdp)
{
	if (parameters.updateSDP_dualityGapThreshold > 0 && duality_gap > parameters.updateSDP_dualityGapThreshold)
	{
		std::cout << "compute_update_sdp : set update_sdp = false \n" << std::flush;
		update_sdp = false;
	}
	else
	{
		std::cout << "compute_update_sdp : set update_sdp = true \n" << std::flush;
		update_sdp = true;
	}
}

void compute_find_zeros(
  const El::BigFloat &duality_gap, const El::BigFloat &primal_objective,
  bool &find_zeros)
{
  find_zeros =  (duality_gap < 0.01);
}  
