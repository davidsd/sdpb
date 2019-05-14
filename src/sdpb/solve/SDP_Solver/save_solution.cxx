#include "../SDP_Solver.hxx"
#include "../../../set_stream_precision.hxx"

#include <boost/filesystem/fstream.hpp>

#include <iomanip>

void SDP_Solver::save_solution(
  const SDP_Solver_Terminate_Reason terminate_reason,
  const std::pair<std::string, Timer> &timer_pair,
  const boost::filesystem::path &out_directory,
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
      out_stream.open(out_directory / "out.txt");
      set_stream_precision(out_stream);
      out_stream << "terminateReason = \"" << terminate_reason << "\";\n"
                 << "primalObjective = " << primal_objective << ";\n"
                 << "dualObjective   = " << dual_objective << ";\n"
                 << "dualityGap      = " << duality_gap << ";\n"
                 << "primalError     = " << primal_error() << ";\n"
                 << "dualError       = " << dual_error << ";\n"
                 << std::setw(16) << std::left << timer_pair.first << "= "
                 << timer_pair.second.elapsed_seconds() << ";\n"
                 << "y = {";
    }
  // y is duplicated among cores, so only need to print out copy on
  // the root node.
  if(!y.blocks.empty())
    {
      El::Print(y.blocks.at(0), "", ",\n", out_stream);
    }

  if(El::mpi::Rank() == 0)
    {
      out_stream << "};\n";
    }

  for(size_t block = 0; block != x.blocks.size(); ++block)
    {
      size_t block_index(block_indices.at(block));
      boost::filesystem::ofstream x_stream(
        out_directory / ("x_" + std::to_string(block_index) + ".txt"));
      El::Print(x.blocks[block], "", "\n", x_stream);
      x_stream << "\n";
    }
}
