#include "../dynamical_solve.hxx"
#include "../set_stream_precision.hxx"
#include "../write_distmatrix.hxx"

#include <boost/filesystem/fstream.hpp>

#include <iomanip>

namespace
{
  void write_psd_block(const boost::filesystem::path &outfile,
                       const El::DistMatrix<El::BigFloat> &block)
  {
    boost::filesystem::ofstream stream;
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
        if(!stream.good())
          {
            throw std::runtime_error("Error when writing to: "
                                     + outfile.string());
          }
      }
  }
}

void save_solution(const Dynamical_Solver &solver,
                   const Dynamical_Solver_Terminate_Reason terminate_reason,
                   const std::pair<std::string, Timer> &timer_pair,
                   const boost::filesystem::path &out_directory,
                   const Write_Solution &write_solution,
                   const std::vector<size_t> &block_indices,
                   const Verbosity &verbosity,
                   const El::Matrix<El::BigFloat> &extParamStep)
{
  // Internally, El::Print() sync's everything to the root core and
  // outputs it from there.  So do not actually open the file on
  // anything but the root node.

  boost::filesystem::ofstream out_stream;
  if(El::mpi::Rank() == 0)
    {
      if(verbosity >= Verbosity::regular)
        {
          std::cout << "Saving solution to      : " << out_directory << '\n';
        }
      const boost::filesystem::path output_path(out_directory / "out.txt");
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
		         << "dualStepSize    = " << solver.d_step << ";\n"
                 << "primalStepSize  = " << solver.p_step << ";\n"
		         << "BFGSHessianUpdated = " << solver.hess_BFGS_updateQ << ";\n"
		         << "NavigatorValue = " << solver.lag_shifted << ";\n"
		         << "findMinimumQ = " << solver.findMinimumQ << ";\n"
		         << "beta = " << solver.final_beta << ";\n"
	         	 << "mulogdetX = " << solver.mulogdetX << ";\n"
                 //<< "totalIterations = " << solver.total_iteration << ";\n"
                 << std::setw(16) << std::left << timer_pair.first << "= "
                 << timer_pair.second.elapsed_seconds() << ";\n";
      if(!out_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + output_path.string());
        }
    }

  boost::filesystem::ofstream extParamStep_stream;
  if(El::mpi::Rank() == 0)
    {
      const boost::filesystem::path extParamStep_path(out_directory / "externalParamStep.txt");
      extParamStep_stream.open(extParamStep_path);
      El::Print(extParamStep,
                std::to_string(extParamStep.Height()) + " "
                + std::to_string(extParamStep.Width()),
               "\n", extParamStep_stream);
      if(!extParamStep_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + extParamStep_path.string());
        } 
    }

  boost::filesystem::ofstream gradient_stream;
  if(El::mpi::Rank() == 0)
    {
      const boost::filesystem::path gradient_path(out_directory / "gradient.txt");
      gradient_stream.open(gradient_path);
      El::Print(solver.grad_BFGS,
                std::to_string(solver.grad_BFGS.Height()) + " "
                + std::to_string(solver.grad_BFGS.Width()),
               "\n", gradient_stream);
      if(!gradient_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + gradient_path.string());
        } 
    }

  boost::filesystem::ofstream gradient_withlog_stream;
  if (El::mpi::Rank() == 0)
  {
	  const boost::filesystem::path gradient_path(out_directory / "gradient_withlog.txt");
	  gradient_withlog_stream.open(gradient_path);
	  El::Print(solver.grad_withlog,
		  std::to_string(solver.grad_withlog.Height()) + " "
		  + std::to_string(solver.grad_withlog.Width()),
		  "\n", gradient_withlog_stream);
	  if (!gradient_withlog_stream.good())
	  {
		  throw std::runtime_error("Error when writing to: "
			  + gradient_path.string());
	  }
  }

  boost::filesystem::ofstream gradient_withoutlog_stream; 
  if (El::mpi::Rank() == 0)
  {
	  const boost::filesystem::path gradient_path(out_directory / "gradient_withoutlog.txt");
	  gradient_withoutlog_stream.open(gradient_path);
	  El::Print(solver.grad_withoutlog,
		  std::to_string(solver.grad_withoutlog.Height()) + " "
		  + std::to_string(solver.grad_withoutlog.Width()),
		  "\n", gradient_withoutlog_stream);
	  if (!gradient_withoutlog_stream.good())
	  {
		  throw std::runtime_error("Error when writing to: "
			  + gradient_path.string());
	  }
  }

  boost::filesystem::ofstream gradient_mixed_stream;
  if (El::mpi::Rank() == 0)
  {
	  const boost::filesystem::path gradient_path(out_directory / "gradient_mixed.txt");
	  gradient_mixed_stream.open(gradient_path);
	  El::Print(solver.grad_mixed,
		  std::to_string(solver.grad_mixed.Height()) + " "
		  + std::to_string(solver.grad_mixed.Width()),
		  "\n", gradient_mixed_stream);
	  if (!gradient_mixed_stream.good())
	  {
		  throw std::runtime_error("Error when writing to: "
			  + gradient_path.string());
	  }
  }

  boost::filesystem::ofstream hess_BFGS_stream;
  if(El::mpi::Rank() == 0)
    {
      const boost::filesystem::path hess_BFGS_path(out_directory / "hessBFGS.txt");
      hess_BFGS_stream.open(hess_BFGS_path);
      El::Print(solver.hess_BFGS,
                std::to_string(solver.hess_BFGS.Height()*solver.hess_BFGS.Width()) + " "
                + std::to_string(1),
               "\n", hess_BFGS_stream);
      if(!hess_BFGS_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + hess_BFGS_path.string());
        } 
    }

  boost::filesystem::ofstream hess_Exact_stream;
  if(El::mpi::Rank() == 0)
    {
      const boost::filesystem::path hess_Exact_path(out_directory / "hessExact.txt");
      hess_Exact_stream.open(hess_Exact_path);
      El::Print(solver.hess_Exact,
                std::to_string(solver.hess_Exact.Height()*solver.hess_Exact.Width()) + " "
                + std::to_string(1),
               "\n", hess_Exact_stream);
      if(!hess_Exact_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + hess_Exact_path.string());
        }
    }

 
  boost::filesystem::ofstream iterations_stream;
  if(El::mpi::Rank() == 0)
    {
      const boost::filesystem::path iterations_path(out_directory / "iterations.txt");
      iterations_stream.open(iterations_path);
      iterations_stream << solver.total_iteration << '\n';
      if(!iterations_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + iterations_path.string());
        }
    }

  //boost::filesystem::ofstream mu_direction_stream;
  //if(El::mpi::Rank() == 0)
  //  {
  //    const boost::filesystem::path mu_direction_path(out_directory / "mu_direction.txt");
  //    mu_direction_stream.open(mu_direction_path);
  //    mu_direction_stream << solver.mu_direction_mode << '\n';
  //    if(!mu_direction_stream.good())
  //      {
  //        throw std::runtime_error("Error when writing to: "
  //                                 + mu_direction_path.string());
  //      }
  //  }

  // y is duplicated among cores, so only need to print out copy on
  // the root node.
  if(write_solution.vector_y && !solver.y.blocks.empty())
    {
      const boost::filesystem::path y_path(out_directory / "y.txt");
      boost::filesystem::ofstream y_stream;
      if(El::mpi::Rank() == 0)
        {
          y_stream.open(y_path);
        }
      El::Print(solver.y.blocks.at(0),
                std::to_string(solver.y.blocks.at(0).Height()) + " "
                  + std::to_string(solver.y.blocks.at(0).Width()),
                "\n", y_stream);
      if(El::mpi::Rank() == 0)
        {
          y_stream << "\n";
          if(!y_stream.good())
            {
              throw std::runtime_error("Error when writing to: "
                                       + y_path.string());
            }
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
