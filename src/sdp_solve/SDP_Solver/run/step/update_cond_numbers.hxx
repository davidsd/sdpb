#pragma once

#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/assert.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdpb_util/cholesky_condition_number.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <El.hpp>

// Calculate and update on rank 0:
// Q_cond_number - condition number Q_cond_number for Cholesky decomposed Q
// max_block_cond_number - max condition number among all blocks of schur_complement_cholesky, X_cholesky,Y_cholesky
// max_block_cond_number_name - name of the corresponding block
inline void update_cond_numbers(
  const El::DistMatrix<El::BigFloat> &Q, const Block_Info &block_info,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Paired_Block_Diagonal_Matrix &X_cholesky,
  const Paired_Block_Diagonal_Matrix &Y_cholesky, Timers &timers,
  // Output variables:
  El::BigFloat &Q_cond_number, El::BigFloat &max_block_cond_number,
  std::string &max_block_cond_number_name)
{
  Scoped_Timer cond_number_timer(timers, "condition_numbers");

  Q_cond_number = cholesky_condition_number(Q);

  max_block_cond_number = 0;
  max_block_cond_number_name = "";

  // schur_complement_cholesky

  const size_t num_blocks = block_info.dimensions.size();
  auto get_schur_block_index = [&](const size_t local_block_index) {
    return block_info.block_indices.at(local_block_index);
  };
  auto [schur_complement_index, schur_complement_cholesky_cond]
    = max_block_cholesky_condition_number(schur_complement_cholesky.blocks,
                                          num_blocks, get_schur_block_index);

  // X_cholesky, Y_cholesky

  const auto num_X_blocks = 2 * num_blocks;
  std::function get_X_index = [&block_info](size_t local_X_index) {
    const size_t parity = local_X_index % 2;
    return 2 * block_info.block_indices.at(local_X_index / 2) + parity;
  };
  auto [X_index, X_cholesky_cond] = max_block_cholesky_condition_number(
    X_cholesky.blocks, num_X_blocks, get_X_index);
  auto [Y_index, Y_cholesky_cond] = max_block_cholesky_condition_number(
    Y_cholesky.blocks, num_X_blocks, get_X_index);

  // Synchronize all at rank 0

  if(El::mpi::Rank() == 0)
    {
      max_block_cond_number = schur_complement_cholesky_cond;
      max_block_cond_number_name = "schur_complement_cholesky.block_"
                                   + std::to_string(schur_complement_index);

      if(max_block_cond_number < X_cholesky_cond)
        {
          max_block_cond_number = X_cholesky_cond;
          const auto block_index = X_index / 2;
          const auto parity = X_index % 2;
          max_block_cond_number_name
            = El::BuildString("X_cholesky.block_", block_index, "_", parity);
        }
      if(max_block_cond_number < Y_cholesky_cond)
        {
          max_block_cond_number = Y_cholesky_cond;
          const auto block_index = Y_index / 2;
          const auto parity = Y_index % 2;
          max_block_cond_number_name
            = El::BuildString("X_cholesky.block_", block_index, "_", parity);
        }
    }
}
