#include "dynamical_solve/dynamical_solve.hxx"
#include "pmp/max_normalization_index.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"
#include "sdpb_util/write_distmatrix.hxx"

#include <fstream>

namespace fs = std::filesystem;

void print_matrix(const El::Matrix<El::BigFloat> &matrix);

namespace
{
  void write_psd_block(const std::filesystem::path &outfile,
                       const El::DistMatrix<El::BigFloat> &block)
  {
    std::ofstream stream;
    if(block.DistRank() == block.Root())
      {
        stream.open(outfile);
      }
    El::Print(block,
              std::to_string(block.Height()) + " "
                + std::to_string(block.Width()),
              "\n", stream);
    if(block.DistRank() == block.Root())
      {
        stream << "\n";
        ASSERT(stream.good(), "Error when writing to: ", outfile);
      }
  }
}

void save_solution(
  const Dynamical_Solver &solver,
  const Dynamical_Solver_Terminate_Reason &terminate_reason,
  const int64_t &solver_runtime, const fs::path &out_directory,
  const Write_Solution &write_solution,
  const std::vector<size_t> &block_indices,
  const std::optional<std::vector<El::BigFloat>> &normalization,
  const Verbosity &verbosity, const El::Matrix<El::BigFloat> &extParamStep)
{
  // Internally, El::Print() sync's everything to the root core and
  // outputs it from there.  So do not actually open the file on
  // anything but the root node.

  std::ofstream out_stream;
  if(El::mpi::Rank() == 0)
    {
      if(verbosity >= Verbosity::regular)
        {
          std::cout << "Saving solution to      : " << out_directory << '\n';
        }
      fs::create_directories(out_directory);
      const fs::path output_path(out_directory / "out.txt");
      out_stream.open(output_path);
      set_stream_precision(out_stream);
      //if (update_sdp && terminate_reason == Dynamical_Solver_Terminate_Reason::MaxIterationsExceeded)
      //  {
      //    out_stream << "exit to update SDPs" << ";\n";
      //  }
      //else
      //  {
      //    out_stream << "terminateReason = \"" << terminate_reason << "\";\n";
      //  }
      out_stream << "terminateReason = \"" << terminate_reason << "\";\n"
                 << "primalObjective = " << solver.primal_objective << ";\n"
                 << "dualObjective   = " << solver.dual_objective << ";\n"
                 << "dualityGap      = " << solver.duality_gap << ";\n"
                 << "primalError     = " << solver.primal_error() << ";\n"
                 << "dualError       = " << solver.dual_error << ";\n"
                 << "Solver runtime  = " << solver_runtime << ";\n";
      ASSERT(out_stream.good(), "Error when writing to: ", output_path);
    }

  std::ofstream out2_stream;
  if(El::mpi::Rank() == 0)
    {
      if(verbosity >= Verbosity::regular)
        {
          std::cout << "Saving solution to      : " << out_directory << '\n';
        }
      const std::filesystem::path output2_path(out_directory
                                               / "skydiving_out.txt");
      out2_stream.open(output2_path);
      set_stream_precision(out2_stream);
      //if (update_sdp && terminate_reason == Dynamical_Solver_Terminate_Reason::MaxIterationsExceeded)
      //  {
      //    out_stream << "exit to update SDPs" << ";\n";
      //  }
      //else
      //  {
      //    out_stream << "terminateReason = \"" << terminate_reason << "\";\n";
      //  }
      out2_stream << "dualStepSize       = " << solver.d_step << ";\n"
                  << "primalStepSize     = " << solver.p_step << ";\n"
                  << "BFGSHessianUpdated = " << solver.hess_BFGS_updateQ
                  << ";\n"
                  << "NavigatorValue     = " << solver.lag_shifted << ";\n"
                  << "findMinimumQ       = " << solver.findMinimumQ << ";\n"
                  << "beta               = " << solver.final_beta << ";\n"
                  << "mulogdetX          = " << solver.mulogdetX << ";\n"
                  << "climbedQ           = " << (!solver.lowest_mu_Q) << ";\n"
                  << "Solver runtime     = " << solver_runtime << ";\n";
      ASSERT(out2_stream.good(), "Error when writing to: ", output2_path);
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

  if(write_solution.vector_y || write_solution.vector_z)
    {
      const fs::path y_path(out_directory / "y.txt");
      const fs::path z_path(out_directory / "z.txt");
      // y is duplicated among blocks, so only need to print out copy
      // from the first block of rank 0.
      const auto &y_dist = solver.y.blocks.at(0);
      if(y_dist.Root() == 0)
        {
          // Copy from all ranks owning a block to rank zero
          El::DistMatrix<El::BigFloat, El::CIRC, El::CIRC> y_circ(y_dist);
          ASSERT_EQUAL(y_circ.Root(), 0);
          if(El::mpi::Rank() == 0)
            {
              // local matrix
              const El::Matrix<El::BigFloat> &y = y_circ.Matrix();
              ASSERT_EQUAL(y.Height(), y_dist.Height());
              ASSERT_EQUAL(y.Width(), 1);
              if(write_solution.vector_y)
                {
                  std::ofstream y_stream(y_path);
                  auto title = El::BuildString(y.Height(), " ", y.Width());
                  El::Print(y, title, "\n", y_stream);
                  y_stream << "\n";
                  ASSERT(y_stream.good(), "Error when writing to: ", y_path);
                }
              if(write_solution.vector_z)
                {
                  ASSERT(normalization.has_value());
                  El::Matrix<El::BigFloat> z(y.Height() + 1, y.Width());
                  ASSERT_EQUAL(normalization->size(), y.Height() + 1);

                  auto max_index
                    = max_normalization_index(normalization.value());
                  ASSERT(normalization->at(max_index) != El::BigFloat(0));
                  // To construct z, we take y
                  // and insert a new element into it at max_index position.
                  // The value of this element is chosen to satisfy
                  // the normalization condition n.z == 1.
                  // NB: this procedure (in particular, choice of max_index)
                  // is the inverse of conversion
                  // Polynomial_Matrix_Program.objective -> Output_SDP.dual_objective_b
                  // in Outout_SDP constructor
                  for(int i = 0; i < max_index; ++i)
                    z(i, 0) = y(i, 0);
                  z(max_index, 0) = 0; // to be updated below
                  for(int i = max_index; i < y.Height(); ++i)
                    z(i + 1, 0) = y(i, 0);

                  // Calculate n.z
                  // we convert sdt::vector to El::Matrix to use it in El::Dot
                  El::Matrix<El::BigFloat> normalization_view;
                  normalization_view.LockedAttach(normalization->size(), 1,
                                                  normalization->data(),
                                                  normalization->size());
                  // n.z - n[max_index]*z[max_index]
                  auto nz = El::Dot(normalization_view, z);
                  z(max_index, 0) = (1 - nz) / normalization->at(max_index);
                  // now n.z == 1

                  // Write to file
                  std::ofstream z_stream(z_path);
                  auto title = El::BuildString(z.Height(), " ", z.Width());
                  El::Print(z, title, "\n", z_stream);
                  z_stream << "\n";
                  ASSERT(z_stream.good(), "Error when writing to: ", z_path);
                }
            }
        }
      if(El::mpi::Rank() == 0)
        {
          // Sanity check.
          // Assertion will fail if something went completely wrong
          // and (y_dist.Root() == 0) is false for all blocks
          if(write_solution.vector_y)
            ASSERT(fs::exists(y_path), y_path);
          if(write_solution.vector_z)
            ASSERT(fs::exists(z_path), z_path);
        }
    }

  for(size_t block = 0; block != solver.x.blocks.size(); ++block)
    {
      size_t block_index(block_indices.at(block));
      if(write_solution.vector_x)
        {
          write_distmatrix(solver.x.blocks.at(block),
                           out_directory
                             / ("x_" + std::to_string(block_index) + ".txt"));
        }
      for(size_t psd_block(0); psd_block < 2; ++psd_block)
        {
          std::string suffix(std::to_string(2 * block_index + psd_block)
                             + ".txt");

          if(write_solution.matrix_X
             && solver.X.blocks.at(2 * block + psd_block).Height() != 0)
            {
              write_psd_block(out_directory / ("X_matrix_" + suffix),
                              solver.X.blocks.at(2 * block + psd_block));
            }
          if(write_solution.matrix_Y
             && solver.Y.blocks.at(2 * block + psd_block).Height() != 0)
            {
              write_psd_block(out_directory / ("Y_matrix_" + suffix),
                              solver.Y.blocks.at(2 * block + psd_block));
            }
        }
    }
}
