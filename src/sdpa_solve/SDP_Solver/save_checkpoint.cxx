#include "../SDP_Solver.hxx"
#include "sdp_solve/SDP_Solver/save_binary_checkpoint.hxx"
#include "sdpb_util/assert.hxx"

#include <filesystem>
#include <boost/property_tree/json_parser.hpp>

namespace fs = std::filesystem;

// We use binary checkpointing because writing text does not write all
// of the necessary digits.  The GMP library sets it to one less than
// required for round-tripping.

namespace Sdpb::Sdpa
{
  void write_binary_checkpoint_data(std::ofstream &checkpoint_stream,
                                    const SDP_Solver &solver)
  {
    // TODO: Write and read precision, num of mpi procs, and procs_per_node.
    write_local_matrix(solver.x, checkpoint_stream);
    write_local_blocks(solver.X, checkpoint_stream);
    write_local_blocks(solver.Y, checkpoint_stream);
  }

  void SDP_Solver::save_checkpoint(
    const fs::path &checkpoint_directory, const Verbosity &verbosity,
    const boost::property_tree::ptree &parameter_properties)
  {
    save_binary_checkpoint(checkpoint_directory, verbosity,
                           parameter_properties, *this);
  }
}
