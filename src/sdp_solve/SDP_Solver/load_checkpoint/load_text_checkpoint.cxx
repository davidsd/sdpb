#include "sdp_solve/SDP_Solver.hxx"
#include "sdp_solve/read_text_block.hxx"

namespace fs = std::filesystem;

bool load_text_checkpoint(const fs::path &checkpoint_directory,
                          const std::vector<size_t> &block_indices,
                          const Verbosity &verbosity, SDP_Solver &solver)
{
  if(!exists(checkpoint_directory / "x_0.txt"))
    {
      return false;
    }

  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << "Loading text checkpoint from : " << checkpoint_directory
                << '\n';
    }

  for(size_t block = 0; block != block_indices.size(); ++block)
    {
      size_t block_index(block_indices.at(block));
      read_text_block(solver.x.blocks.at(block), checkpoint_directory, "x_",
                      block_index);
      read_text_block(solver.y.blocks.at(block),
                      checkpoint_directory / "y.txt");

      for(size_t psd_block(0); psd_block < 2; ++psd_block)
        {
          // Constant constraints have empty odd parity blocks, so we do not
          // need to load them.
          if(solver.X.blocks.at(2 * block + psd_block).Height() != 0)
            {
              const size_t psd_index(2 * block_index + psd_block);
              read_text_block(solver.X.blocks.at(2 * block + psd_block),
                              checkpoint_directory, "X_matrix_", psd_index);
              read_text_block(solver.Y.blocks.at(2 * block + psd_block),
                              checkpoint_directory, "Y_matrix_", psd_index);
            }
        }
    }
  return true;
}
