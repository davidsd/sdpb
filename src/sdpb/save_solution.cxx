#include "../sdp_solve.hxx"
#include "../set_stream_precision.hxx"
#include "../write_distmatrix.hxx"

namespace fs = std::filesystem;

namespace
{
  void write_psd_block(const fs::path &outfile,
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
        if(!stream.good())
          {
            throw std::runtime_error("Error when writing to: "
                                     + outfile.string());
          }
      }
  }
}

void save_solution(const SDP_Solver &solver,
                   const SDP_Solver_Terminate_Reason &terminate_reason,
                   const int64_t &solver_runtime,
                   const fs::path &out_directory,
                   const Write_Solution &write_solution,
                   const std::vector<size_t> &block_indices,
                   const Verbosity &verbosity)
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
      out_stream << "terminateReason = \"" << terminate_reason << "\";\n"
                 << "primalObjective = " << solver.primal_objective << ";\n"
                 << "dualObjective   = " << solver.dual_objective << ";\n"
                 << "dualityGap      = " << solver.duality_gap << ";\n"
                 << "primalError     = " << solver.primal_error() << ";\n"
                 << "dualError       = " << solver.dual_error << ";\n"
                 << "Solver runtime  = " << solver_runtime << ";\n";
      if(!out_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + output_path.string());
        }
    }
  // y is duplicated among cores, so only need to print out copy on
  // the root node.
  if(write_solution.vector_y && !solver.y.blocks.empty())
    {
      const fs::path y_path(out_directory / "y.txt");
      std::ofstream y_stream;
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
