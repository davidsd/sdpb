#include "save_binary_checkpoint.hxx"
#include "../SDP_Solver.hxx"
#include "sdpb_util/assert.hxx"

#include <filesystem>
#include <boost/property_tree/json_parser.hpp>

namespace fs = std::filesystem;

void write_binary_checkpoint_data(std::ofstream &checkpoint_stream,
                                  const SDP_Solver &solver)
{
  // TODO: Write and read precision, num of mpi procs, and procs_per_node.
  write_local_blocks(solver.x, checkpoint_stream);
  write_local_blocks(solver.X, checkpoint_stream);
  write_local_blocks(solver.y, checkpoint_stream);
  write_local_blocks(solver.Y, checkpoint_stream);
}

void SDP_Solver::save_checkpoint(
  const fs::path &checkpoint_directory, const Verbosity &verbosity,
  const boost::property_tree::ptree &parameter_properties)
{
  save_binary_checkpoint(checkpoint_directory, verbosity, parameter_properties,
                         *this);
}
