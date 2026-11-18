#pragma once

#include "sdp_solve/SDP_Solver.hxx"
#include "sdp_solve/read_text_block.hxx"
#include "sdpb_util/assert.hxx"

#include <string>


inline bool load_text_checkpoint(const std::filesystem::path &checkpoint_directory,
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

  // X and Y are read for all blocks or for none (a checkpoint with x and y
  // only, for sdpb --reconstructXYmu); some of them is an error.
  int has_XY = 0, lacks_XY = 0;
  for(size_t block = 0; block != block_indices.size(); ++block)
    for(size_t psd_block(0); psd_block < 2; ++psd_block)
      {
        const size_t psd_index = 2 * block_indices.at(block) + psd_block;
        if(solver.X.get_block(block, psd_block).Height() == 0)
          continue;
        for(const std::string prefix : {"X_matrix_", "Y_matrix_"})
          {
            if(exists(checkpoint_directory
                      / (prefix + std::to_string(psd_index) + ".txt")))
              has_XY = 1;
            else
              lacks_XY = 1;
          }
      }
  has_XY = El::mpi::AllReduce(has_XY, El::mpi::MAX, El::mpi::COMM_WORLD);
  lacks_XY = El::mpi::AllReduce(lacks_XY, El::mpi::MAX, El::mpi::COMM_WORLD);
  if(has_XY && lacks_XY)
    RUNTIME_ERROR("Text checkpoint ", checkpoint_directory,
                  " has X_matrix/Y_matrix files for some blocks but not "
                  "for others.");
  solver.loaded_without_XY = !has_XY;

  for(size_t block = 0; block != block_indices.size(); ++block)
    {
      size_t block_index(block_indices.at(block));
      read_text_block(solver.x.blocks.at(block), checkpoint_directory, "x_",
                      block_index);
      read_text_block(solver.y.blocks.at(block),
                      checkpoint_directory / "y.txt");

      for(size_t psd_block(0); psd_block < 2; ++psd_block)
        {
          const size_t psd_index = 2 * block_index + psd_block;
          auto &X_block = solver.X.get_block(block, psd_block);
          auto &Y_block = solver.Y.get_block(block, psd_block);
          // Constant constraints have empty odd parity blocks, so we do not
          // need to load them. A checkpoint with only x and y has no X/Y
          // files: X and Y stay 0, see sdpb --reconstructXYmu.
          if(X_block.Height() != 0 && has_XY)
            {
              read_text_block(X_block, checkpoint_directory, "X_matrix_",
                              psd_index);
              read_text_block(Y_block, checkpoint_directory, "Y_matrix_",
                              psd_index);
            }
        }
    }
  return true;
}
