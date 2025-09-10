#include "load_binary_checkpoint.hxx"
#include "load_text_checkpoint.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdpb_util/assert.hxx"

namespace fs = std::filesystem;

void load_local_binary_data(SDP_Solver &solver,
                            std::ifstream &checkpoint_stream)
{
  read_local_binary_blocks(solver.x, checkpoint_stream);
  read_local_binary_blocks(solver.X, checkpoint_stream);
  read_local_binary_blocks(solver.y, checkpoint_stream);
  read_local_binary_blocks(solver.Y, checkpoint_stream);
}

bool SDP_Solver::load_checkpoint(const fs::path &checkpoint_directory,
                                 const Block_Info &block_info,
                                 const Verbosity &verbosity,
                                 const bool &require_initial_checkpoint)
{
  bool valid_checkpoint(
    load_binary_checkpoint(checkpoint_directory, verbosity, *this)
    || load_text_checkpoint(checkpoint_directory, block_info.block_indices,
                            verbosity, *this));
  if(!valid_checkpoint && require_initial_checkpoint)
    {
      RUNTIME_ERROR("Unable to load checkpoint from directory: ",
                    checkpoint_directory);
    }
  return valid_checkpoint;
}
