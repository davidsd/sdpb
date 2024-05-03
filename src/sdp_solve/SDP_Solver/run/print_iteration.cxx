#include "sdp_solve/SDP_Solver.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <chrono>
#include <iostream>
#include <iomanip>

void print_iteration(
  const std::filesystem::path &iterations_json_path, const int &iteration,
  const El::BigFloat &mu, const El::BigFloat &primal_step_length,
  const El::BigFloat &dual_step_length, const El::BigFloat &beta_corrector,
  const SDP_Solver &sdp_solver,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &iteration_start_time,
  const El::BigFloat &Q_cond_number, const El::BigFloat &max_block_cond_number,
  const std::string &max_block_cond_number_name, const Verbosity &verbosity)
{
  if(El::mpi::Rank() != 0)
    return;

  const double iteration_time_seconds
    = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - iteration_start_time)
        .count()
      / 1000.0;
  const double runtime_seconds
    = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - solver_start_time)
        .count()
      / 1000.0;

  // Print to stdout
  if(verbosity >= Verbosity::regular)
    {
      std::cout << std::left << std::setw(4) << iteration << "  "

                << std::right << std::setw(8)
                << static_cast<size_t>(runtime_seconds) << " "

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

  // Add iteration to out/iterations.json
  // TODO: always or only for verbosity >= regular?
  // Seems that we can write it always, if it isn't too slow.
  std::ofstream os_json(iterations_json_path, std::ios::app);
  if(os_json.good())
    {
      if(iteration != 1)
        os_json << ",";

      os_json << "\n{ \"iteration\":" << iteration;
      os_json << std::setprecision(3) << std::fixed;
      os_json << ", \"total_time\": " << runtime_seconds
              << ", \"iter_time\": " << iteration_time_seconds;
      os_json << std::defaultfloat;
      set_stream_precision(os_json);
      os_json << ", \"mu\": \"" << mu << "\""
              << ", \"P-obj\": \"" << sdp_solver.primal_objective << "\""
              << ", \"D-obj\": \"" << sdp_solver.dual_objective << "\""
              << ", \"gap\": \"" << sdp_solver.duality_gap << "\""
              << ", \"P-err\": \"" << sdp_solver.primal_error_P << "\""
              << ", \"p-err\": \"" << sdp_solver.primal_error_p << "\""
              << ", \"D-err\": \"" << sdp_solver.dual_error << "\""
              << ", \"R-err\": \"" << sdp_solver.R_error << "\""
              << ", \"P-step\": \"" << primal_step_length << "\""
              << ", \"D-step\": \"" << dual_step_length << "\""
              << ", \"beta\": \"" << beta_corrector << "\""
              << ", \"Q_cond_number\": \"" << Q_cond_number << "\""
              << ", \"max_block_cond_number\": \"" << max_block_cond_number
              << "\""
              << ", \"block_name\": \"" << max_block_cond_number_name << "\""
              << " }";
      ASSERT(Q_cond_number >= 1, DEBUG_STRING(Q_cond_number));
      ASSERT(max_block_cond_number >= 1, DEBUG_STRING(max_block_cond_number));
      ASSERT(!max_block_cond_number_name.empty());
    }
  if(verbosity >= Verbosity::trace)
    {
      El::Output("Cholesky condition number of Q: ", Q_cond_number);
      El::Output("Max Cholesky condition number among blocks: ",
                 max_block_cond_number, ", ", max_block_cond_number_name);
    }
}
