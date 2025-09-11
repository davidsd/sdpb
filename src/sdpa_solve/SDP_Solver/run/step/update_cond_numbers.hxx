#pragma once

#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpa_solve/Block_Info.hxx"
#include "sdpb_util/cholesky_condition_number.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <string>

namespace Sdpb::Sdpa
{
  // Calculate and update on rank 0:
  // S_cond_number - condition number S_cond_number for Cholesky decomposed S
  // max_block_cond_number - max condition number among all blocks of schur_complement_cholesky, X_cholesky,Y_cholesky
  // max_block_cond_number_name - name of the corresponding block
  inline void
  update_cond_numbers(const El::DistMatrix<El::BigFloat> &S,
                      const Block_Info &block_info,
                      const Block_Diagonal_Matrix &X_cholesky,
                      const Block_Diagonal_Matrix &Y_cholesky, Timers &timers,
                      // Output variables:
                      El::BigFloat &S_cond_number,
                      El::BigFloat &max_block_cond_number,
                      std::string &max_block_cond_number_name)
  {
    Scoped_Timer cond_number_timer(timers, "condition_numbers");

    S_cond_number = cholesky_condition_number(S);

    max_block_cond_number = 0;
    max_block_cond_number_name = "";

    const size_t num_blocks = block_info.block_dimensions.size();

    auto get_block_index = [&](const size_t local_block_index) {
      return block_info.block_indices.at(local_block_index);
    };

    // X_cholesky, Y_cholesky

    auto [X_index, X_cholesky_cond] = max_block_cholesky_condition_number(
      X_cholesky.blocks, num_blocks, get_block_index);
    auto [Y_index, Y_cholesky_cond] = max_block_cholesky_condition_number(
      Y_cholesky.blocks, num_blocks, get_block_index);

    if(El::mpi::Rank() == 0)
      {
        max_block_cond_number = X_cholesky_cond;
        max_block_cond_number_name
          = El::BuildString("X_cholesky.block_", X_index);
        if(max_block_cond_number < Y_cholesky_cond)
          {
            max_block_cond_number = Y_cholesky_cond;
            max_block_cond_number_name
              = El::BuildString("X_cholesky.block_", Y_index);
          }
      }
  }
}
