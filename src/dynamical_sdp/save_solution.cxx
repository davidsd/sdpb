#include "dynamical_solve/dynamical_solve.hxx"
#include "pmp/max_normalization_index.hxx"
#include "sdpb/save_sdpb_solution.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <fstream>

namespace fs = std::filesystem;

void print_matrix(const El::Matrix<El::BigFloat> &matrix);

void save_solution(
  const Dynamical_Solver &solver,
  const Dynamical_Solver_Terminate_Reason &terminate_reason,
  const int64_t &solver_runtime, const fs::path &out_directory,
  const Write_Solution &write_solution,
  const std::vector<size_t> &block_indices,
  const std::optional<std::vector<El::BigFloat>> &normalization,
  const Verbosity &verbosity, const El::Matrix<El::BigFloat> &extParamStep)
{
  // Regular SDPB output
  save_sdpb_solution(solver, terminate_reason, solver_runtime, out_directory,
                     write_solution, block_indices, normalization, verbosity);

  // Skydiving output

  if(El::mpi::Rank() == 0)
    {
      std::ofstream skydiving_out_stream;
      if(verbosity >= Verbosity::regular)
        {
          std::cout << "Saving solution to      : " << out_directory << '\n';
        }
      const std::filesystem::path skydiving_out_path(out_directory
                                                     / "skydiving_out.txt");
      skydiving_out_stream.open(skydiving_out_path);
      set_stream_precision(skydiving_out_stream);
      //if (update_sdp && terminate_reason == Dynamical_Solver_Terminate_Reason::MaxIterationsExceeded)
      //  {
      //    skydiving_out_stream << "exit to update SDPs" << ";\n";
      //  }
      //else
      //  {
      //    skydiving_out_stream << "terminateReason = \"" << terminate_reason << "\";\n";
      //  }
      skydiving_out_stream
        << "dualStepSize       = " << solver.d_step << ";\n"
        << "primalStepSize     = " << solver.p_step << ";\n"
        << "BFGSHessianUpdated = " << solver.hess_BFGS_updateQ << ";\n"
        << "NavigatorValue     = " << solver.lag_shifted << ";\n"
        << "findMinimumQ       = " << solver.findMinimumQ << ";\n"
        << "beta               = " << solver.final_beta << ";\n"
        << "mulogdetX          = " << solver.mulogdetX << ";\n"
        << "climbedQ           = " << (!solver.lowest_mu_Q) << ";\n"
        << "Solver runtime     = " << solver_runtime << ";\n";
      ASSERT(skydiving_out_stream.good(),
             "Error when writing to: ", skydiving_out_path);
    }

  std::ofstream extParamStep_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path extParamStep_path(out_directory
                                                    / "externalParamStep.txt");
      extParamStep_stream.open(extParamStep_path);
      El::Print(extParamStep,
                std::to_string(extParamStep.Height()) + " "
                  + std::to_string(extParamStep.Width()),
                "\n", extParamStep_stream);
      ASSERT(extParamStep_stream.good(),
             "Error when writing to: ", extParamStep_path);
    }

  std::ofstream gradient_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path gradient_path(out_directory
                                                / "gradient.txt");
      gradient_stream.open(gradient_path);
      El::Print(solver.grad_BFGS,
                std::to_string(solver.grad_BFGS.Height()) + " "
                  + std::to_string(solver.grad_BFGS.Width()),
                "\n", gradient_stream);
      ASSERT(gradient_stream.good(), "Error when writing to: ", gradient_path);
    }

  // SN_V15
  std::ofstream gradient_withlog_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path gradient_path(out_directory
                                                / "gradient_withlog.txt");
      gradient_withlog_stream.open(gradient_path);
      El::Print(solver.grad_withlog,
                std::to_string(solver.grad_withlog.Height()) + " "
                  + std::to_string(solver.grad_withlog.Width()),
                "\n", gradient_withlog_stream);
      ASSERT(gradient_withlog_stream.good(),
             "Error when writing to: ", gradient_path);
    }

  std::ofstream gradient_withoutlog_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path gradient_path(out_directory
                                                / "gradient_withoutlog.txt");
      gradient_withoutlog_stream.open(gradient_path);
      El::Print(solver.grad_withoutlog,
                std::to_string(solver.grad_withoutlog.Height()) + " "
                  + std::to_string(solver.grad_withoutlog.Width()),
                "\n", gradient_withoutlog_stream);
      ASSERT(gradient_withoutlog_stream.good(),
             "Error when writing to: ", gradient_path);
    }

  std::ofstream gradient_mixed_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path gradient_path(out_directory
                                                / "gradient_mixed.txt");
      gradient_mixed_stream.open(gradient_path);
      El::Print(solver.grad_mixed,
                std::to_string(solver.grad_mixed.Height()) + " "
                  + std::to_string(solver.grad_mixed.Width()),
                "\n", gradient_mixed_stream);
      ASSERT(gradient_mixed_stream.good(),
             "Error when writing to: ", gradient_path);
    }

  // SN_V15
  std::ofstream hess_BFGS_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path hess_BFGS_path(out_directory
                                                 / "hessBFGS.txt");
      hess_BFGS_stream.open(hess_BFGS_path);
      El::Print(
        solver.hess_BFGS,
        std::to_string(solver.hess_BFGS.Height() * solver.hess_BFGS.Width())
          + " " + std::to_string(1),
        "\n", hess_BFGS_stream);
      ASSERT(hess_BFGS_stream.good(),
             "Error when writing to: ", hess_BFGS_path);
    }

  std::ofstream hess_pp_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path hess_pp_path(out_directory
                                               / "hessBFGSpp.txt");
      hess_pp_stream.open(hess_pp_path);
      El::Print(solver.hess_BFGS_pp,
                std::to_string(solver.hess_BFGS_pp.Height()
                               * solver.hess_BFGS_pp.Width())
                  + " " + std::to_string(1),
                "\n", hess_pp_stream);
      ASSERT(hess_pp_stream.good(), "Error when writing to: ", hess_pp_path);
    }

  std::ofstream hess_Exact_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path hess_Exact_path(out_directory
                                                  / "hessExact.txt");
      hess_Exact_stream.open(hess_Exact_path);
      El::Print(
        solver.hess_Exact,
        std::to_string(solver.hess_Exact.Height() * solver.hess_Exact.Width())
          + " " + std::to_string(1),
        "\n", hess_Exact_stream);
      ASSERT(hess_Exact_stream.good(),
             "Error when writing to: ", hess_Exact_path);
    }

  std::ofstream iterations_stream;
  if(El::mpi::Rank() == 0)
    {
      const std::filesystem::path iterations_path(out_directory
                                                  / "iterations.txt");
      iterations_stream.open(iterations_path);
      iterations_stream << solver.total_iteration << '\n';
      ASSERT(iterations_stream.good(),
             "Error when writing to: ", iterations_path);
    }

  // std::ofstream mu_direction_stream;
  // if(El::mpi::Rank() == 0)
  //   {
  //     const std::filesystem::path mu_direction_path(out_directory
  //                                                   / "mu_direction.txt");
  //     mu_direction_stream.open(mu_direction_path);
  //     mu_direction_stream << solver.mu_direction_mode << '\n';
  //     ASSERT(mu_direction_stream.good(),
  //            "Error when writing to: ", mu_direction_path);
  //   }
}
