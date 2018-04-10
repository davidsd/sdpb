#include "../Timers.hxx"
#include "../SDP_Solver.hxx"

#include <boost/date_time/posix_time/posix_time.hpp>

#include <iostream>

void SDP_Solver::print_iteration(int iteration, Real mu,
                                 Real primal_step_length,
                                 Real dual_step_length, Real beta_corrector)
{
  boost::posix_time::time_duration td(
    boost::posix_time::microseconds(timers["Solver runtime"].elapsed().wall)
    / 1000);
  std::stringstream ss;
  ss << td;
  gmp_fprintf(stdout,
              "%3d  %s  %-8.1Fe %-+11.2Fe %-+11.2Fe %-9.2Fe  %-+10.2Fe  "
              "%-+10.2Fe  %-8.3Fg %-8.3Fg %-4.2Fg  %d",
              iteration, ss.str().substr(0, 8).c_str(), mu.get_mpf_t(),
              primal_objective.get_mpf_t(), dual_objective.get_mpf_t(),
              duality_gap.get_mpf_t(), primal_error.get_mpf_t(),
              dual_error.get_mpf_t(), primal_step_length.get_mpf_t(),
              dual_step_length.get_mpf_t(), beta_corrector.get_mpf_t(),
              static_cast<int>(sdp.dual_objective.size()));
  std::cout << '\n';
}
