#include "dynamical_solve/Dynamical_Solver.hxx"

bool load_binary_checkpoint(const std::filesystem::path &checkpoint_directory,
                            const Verbosity &verbosity,
                            Dynamical_Solver &solver);

bool load_text_checkpoint(const std::filesystem::path &checkpoint_directory,
                          const std::vector<size_t> &block_indices,
                          const Verbosity &verbosity,
                          Dynamical_Solver &solver);

bool Dynamical_Solver::load_checkpoint(
  const std::filesystem::path &checkpoint_directory,
  const Block_Info &block_info, const Verbosity &verbosity,
  const bool &require_initial_checkpoint)
{
  bool valid_checkpoint(
    load_binary_checkpoint(checkpoint_directory, verbosity, *this)
    || load_text_checkpoint(checkpoint_directory, block_info.block_indices,
                            verbosity, *this));
  if(!valid_checkpoint && require_initial_checkpoint)
    {
      throw std::runtime_error("Unable to load checkpoint from directory: "
                               + checkpoint_directory.string());
    }
  return valid_checkpoint;
}
