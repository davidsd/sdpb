#include "../SDP_Solver.hxx"
#include "../../Timers.hxx"

#include <boost/date_time/posix_time/posix_time.hpp>

#include <iostream>

void SDP_Solver::print_iteration(const int &iteration, El::BigFloat &mu,
                                 const El::BigFloat &primal_step_length,
                                 const El::BigFloat &dual_step_length,
                                 const El::BigFloat &beta_corrector,
                                 const size_t &dual_objective_b_height)
{
  if(El::mpi::Rank() == 0)
    {
      boost::posix_time::time_duration td(
        boost::posix_time::microseconds(
          timers["Solver runtime"].elapsed().wall)
        / 1000);

      std::stringstream ss;
      ss << td;

      std::cout << std::left << std::setw(3) << iteration << "  "
                << ss.str().substr(0, 8) << "  "

                << std::left << std::setw(8) << std::setprecision(2)
                << std::showpoint << static_cast<double>(mu) << " "

                << std::showpos << std::setw(11) << std::setprecision(3)
                << static_cast<double>(primal_objective) << " "

                << std::setw(11) << std::setprecision(3)
                << static_cast<double>(dual_objective) << " "

                << std::noshowpos << std::setw(10) << std::setprecision(3)
                << static_cast<double>(duality_gap) << " "

                << std::showpos << std::setw(11) << std::setprecision(3)
                << static_cast<double>(primal_error) << " "

                << std::setw(11) << std::setprecision(3)
                << static_cast<double>(dual_error) << " "

                << std::noshowpos << std::setw(8) << std::setprecision(3)
                << static_cast<double>(primal_step_length) << " "

                << std::setw(8) << std::setprecision(3)
                << static_cast<double>(dual_step_length) << " "

                << std::setw(4) << std::setprecision(3)
                << static_cast<double>(beta_corrector) << "  "

                << dual_objective_b_height << "\n";
    }
}
