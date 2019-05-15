#include "../../SDP_Solver.hxx"

bool load_binary_checkpoint(const boost::filesystem::path &checkpoint_directory,
                            const Verbosity &verbosity, SDP_Solver &solver);

bool load_text_checkpoint(const boost::filesystem::path &checkpoint_directory,
                          const std::vector<size_t> &block_indices,
                          const Verbosity &verbosity, SDP_Solver &solver);

bool SDP_Solver::load_checkpoint(
  const boost::filesystem::path &checkpoint_directory,
  const Block_Info &block_info, const Verbosity &verbosity)
{
  return load_binary_checkpoint(checkpoint_directory, verbosity, *this)
         || load_text_checkpoint(checkpoint_directory,
                                 block_info.block_indices, verbosity, *this);
}
