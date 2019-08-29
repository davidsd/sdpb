#include "../SDP_Solver.hxx"
#include "../../../set_stream_precision.hxx"

#include <boost/filesystem/fstream.hpp>

#include <iomanip>

namespace
{
  void write_psd_block(const boost::filesystem::path &outfile,
                       const El::DistMatrix<El::BigFloat> &block)
  {
    boost::filesystem::ofstream stream(outfile);
    El::Print(block,
              std::to_string(block.Height()) + " "
                + std::to_string(block.Width()),
              "\n", stream);
    stream << "\n";
    if(!stream.good())
      {
        throw std::runtime_error("Error when writing to: " + outfile.string());
      }
  }
}

void SDP_Solver::save_solution(
  const SDP_Solver_Terminate_Reason terminate_reason,
  const std::pair<std::string, Timer> &timer_pair,
  const boost::filesystem::path &out_directory, const bool &write_matrices,
  const std::vector<size_t> &block_indices, const Verbosity &verbosity) const
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
      out_stream << "terminateReason = \"" << terminate_reason << "\";\n"
                 << "primalObjective = " << primal_objective << ";\n"
                 << "dualObjective   = " << dual_objective << ";\n"
                 << "dualityGap      = " << duality_gap << ";\n"
                 << "primalError     = " << primal_error() << ";\n"
                 << "dualError       = " << dual_error << ";\n"
                 << std::setw(16) << std::left << timer_pair.first << "= "
                 << timer_pair.second.elapsed_seconds() << ";\n";
      if(!out_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + output_path.string());
        }
    }
  // y is duplicated among cores, so only need to print out copy on
  // the root node.
  if(!y.blocks.empty())
    {
      const boost::filesystem::path y_path(out_directory / "y.txt");
      boost::filesystem::ofstream y_stream(y_path);
      El::Print(y.blocks.at(0),
                std::to_string(y.blocks.at(0).Height()) + " "
                  + std::to_string(y.blocks.at(0).Width()),
                "\n", y_stream);
      if(!y_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + y_path.string());
        }
    }

  for(size_t block = 0; block != x.blocks.size(); ++block)
    {
      size_t block_index(block_indices.at(block));
      const boost::filesystem::path x_path(
        out_directory / ("x_" + std::to_string(block_index) + ".txt"));
      boost::filesystem::ofstream x_stream(x_path);
      El::Print(x.blocks.at(block),
                std::to_string(x.blocks.at(block).Height()) + " "
                  + std::to_string(x.blocks.at(block).Width()),
                "\n", x_stream);
      if(!x_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + x_path.string());
        }

      if(write_matrices)
        {
          for(size_t psd_block(0); psd_block < 2; ++psd_block)
            {
              std::string suffix(std::to_string(2 * block_index + psd_block)
                                 + ".txt");

              write_psd_block(out_directory / ("X_matrix_" + suffix),
                              X.blocks.at(2 * block + psd_block));
              write_psd_block(out_directory / ("Y_matrix_" + suffix),
                              Y.blocks.at(2 * block + psd_block));
            }
        }
    }
}
