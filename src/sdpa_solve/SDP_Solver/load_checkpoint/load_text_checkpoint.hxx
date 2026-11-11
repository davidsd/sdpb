#pragma once

#include "sdpa_solve/SDP_Solver.hxx"
#include "sdp_solve/read_text_block.hxx"

namespace Sdpb::Sdpa
{
  inline bool load_text_checkpoint(const std::filesystem::path &checkpoint_directory,
                            const std::vector<size_t> &block_indices,
                            const Verbosity &verbosity, SDP_Solver &solver)
  {
    if(!exists(checkpoint_directory / "x.txt"))
      {
        return false;
      }

    if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
      {
        std::cout << "Loading text checkpoint from : " << checkpoint_directory
                  << '\n';
      }

    read_text_block(solver.x, checkpoint_directory / "x.txt");
    for(size_t block = 0; block != block_indices.size(); ++block)
      {
        const auto block_index = block_indices.at(block);
        if(solver.X.blocks.at(block).Height() != 0)
          {
            read_text_block(solver.X.blocks.at(block), checkpoint_directory,
                            "X_matrix_", block_index);
            read_text_block(solver.Y.blocks.at(block), checkpoint_directory,
                            "Y_matrix_", block_index);
          }
      }
    return true;
  }
}