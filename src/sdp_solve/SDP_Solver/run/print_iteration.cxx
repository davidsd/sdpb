#include "../../SDP_Solver.hxx"

#include <chrono>
#include <iostream>
#include <iomanip>

void print_iteration(
  const int &iteration, const El::BigFloat &mu,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const El::BigFloat &beta_corrector, const SDP_Solver &sdp_solver,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
  &solver_start_time,
                     const Verbosity &verbosity)
{
  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << std::left << std::setw(4) << iteration << "  "

                << std::right << std::setw(8)
                << std::chrono::duration_cast<std::chrono::seconds>(
                     std::chrono::high_resolution_clock::now()
                     - solver_start_time)
                     .count()
                << " "

                << std::left << std::setw(8) << std::setprecision(2)
                << std::showpoint << static_cast<double>(mu) << " "

                << std::showpos << std::setw(11) << std::setprecision(3)
                << static_cast<double>(sdp_solver.primal_objective) << " "

                << std::setw(11) << std::setprecision(3)
                << static_cast<double>(sdp_solver.dual_objective) << " "

                << std::noshowpos << std::setw(10) << std::setprecision(3)
                << static_cast<double>(sdp_solver.duality_gap) << " "

                << std::showpos << std::setw(11) << std::setprecision(3)
                << static_cast<double>(sdp_solver.primal_error_P) << " "
                << std::showpos << std::setw(11) << std::setprecision(3)
                << static_cast<double>(sdp_solver.primal_error_p) << " "

                << std::setw(11) << std::setprecision(3)
                << static_cast<double>(sdp_solver.dual_error) << " "

                << std::noshowpos << std::setw(8) << std::setprecision(3)
                << static_cast<double>(primal_step_length) << " "

                << std::setw(8) << std::setprecision(3)
                << static_cast<double>(dual_step_length) << " "

                << std::setw(4) << std::setprecision(3)
                << static_cast<double>(beta_corrector) << "\n"

                << std::flush;
    }
}
