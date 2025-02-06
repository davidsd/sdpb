#include "sdp_solve/SDP_Solver.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"
#include "step/Corrector_Iteration.hxx"

#include <chrono>
#include <iostream>
#include <iomanip>

void print_iteration(
  const std::filesystem::path &iterations_json_path, const int &iteration,
  const El::BigFloat &mu, const El::BigFloat &primal_step_length,
  const El::BigFloat &dual_step_length, const El::BigFloat &beta_corrector,
  const std::vector<Corrector_Iteration> &corrector_iterations,
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
      os_json << "\n{";

      bool need_comma = false;
      const auto add_value
        = [&os_json, &need_comma](const auto &name, const auto &value) {
            if(need_comma)
              os_json << ",";
            os_json << "\"" << name << "\":" << value;
            need_comma = true;
          };

      const auto add_string
        = [&os_json, &need_comma](const auto &name, const auto &value) {
            if(need_comma)
              os_json << ",";
            os_json << "\"" << name << "\": \"" << value << "\"";
            need_comma = true;
          };

      add_value("iteration", iteration);
      os_json << std::setprecision(3) << std::fixed;

      add_value("total_time", runtime_seconds);
      add_value("iter_time", iteration_time_seconds);
      add_value("num_corrector_iterations", corrector_iterations.size());

      os_json << std::defaultfloat;
      set_stream_precision(os_json);

      add_string("mu", mu);
      add_string("P-obj", sdp_solver.primal_objective);
      add_string("D-obj", sdp_solver.dual_objective);
      add_string("gap", sdp_solver.duality_gap);
      add_string("P-err", sdp_solver.primal_error_P);
      add_string("p-err", sdp_solver.primal_error_p);
      add_string("D-err", sdp_solver.dual_error);
      add_string("R-err", sdp_solver.R_error);
      add_string("P-step", primal_step_length);
      add_string("D-step", dual_step_length);
      add_string("beta", beta_corrector);
      add_string("Q_cond_number", Q_cond_number);
      add_string("max_block_cond_number", max_block_cond_number);
      add_string("block_name", max_block_cond_number_name);
      os_json << ",\"corrector_iterations\": [";
      for(size_t i = 0; i < corrector_iterations.size(); i++)
        {
          if(i != 0)
            os_json << ",";
          os_json << "{";
          {
            const auto &iter = corrector_iterations.at(i);
            need_comma = false;
            add_string("P-step", iter.primal_step_length);
            add_string("D-step", iter.dual_step_length);
            add_string("max P-step", iter.max_primal_step_length);
            add_string("max D-step", iter.max_dual_step_length);
            add_string("mu", iter.mu);
            add_string("R-err", iter.R_error);
            add_string("R_mean_abs", iter.R_mean_abs);
            add_string("dXdY_mean", iter.dXdY_mean);
            add_string("-dlog(mu)/dt_full", iter.log_mu_speed_full);
            add_string("-dlog(mu)/dt_corr", iter.log_mu_speed_corrector);
            add_value("is_canceled", iter.is_canceled ? "true" : "false");
          }
          os_json << "}";
        }
      os_json << "]";
      os_json << " }";
      ASSERT(Q_cond_number >= 1, DEBUG_STRING(Q_cond_number));
      ASSERT(max_block_cond_number >= 1, DEBUG_STRING(max_block_cond_number));
      ASSERT(!max_block_cond_number_name.empty());

#undef ADD_VALUE
#undef ADD_QUOTED_NO_COMMA
#undef ADD_QUOTED
    }
  if(verbosity >= Verbosity::trace)
    {
      El::Output("Cholesky condition number of Q: ", Q_cond_number);
      El::Output("Max Cholesky condition number among blocks: ",
                 max_block_cond_number, ", ", max_block_cond_number_name);
    }
}
