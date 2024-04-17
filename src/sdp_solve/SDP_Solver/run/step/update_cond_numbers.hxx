#pragma once

#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/assert.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdpb_util/cholesky_condition_number.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <El.hpp>

// Calculate and update on rank 0:
// Q_cond_number - condition number Q_cond_number for Cholesky decomposed Q
// max_block_cond_number - max condition number among all blocks of schur_complement_cholesky, X_cholesky,Y_cholesky
// max_block_cond_number_name - name of the corresponding block
inline void
update_cond_numbers(const El::DistMatrix<El::BigFloat> &Q,
                    const Block_Info &block_info,
                    const Block_Diagonal_Matrix &schur_complement_cholesky,
                    const Block_Diagonal_Matrix &X_cholesky,
                    const Block_Diagonal_Matrix &Y_cholesky, Timers &timers,
                    // Output variables:
                    El::BigFloat &Q_cond_number,
                    El::BigFloat &max_block_cond_number,
                    std::string &max_block_cond_number_name)
{
  Scoped_Timer cond_number_timer(timers, "condition_numbers");

  Q_cond_number = cholesky_condition_number(Q);

  max_block_cond_number = 0;
  max_block_cond_number_name = "";

  // schur_complement_cholesky

  const size_t num_blocks = block_info.dimensions.size();
  std::vector<El::BigFloat> schur_complement_cholesky_cond(num_blocks);
  for(size_t local_block_index = 0;
      local_block_index < schur_complement_cholesky.blocks.size();
      ++local_block_index)
    {
      const auto &cholesky_block
        = schur_complement_cholesky.blocks.at(local_block_index);
      const auto global_block_index
        = block_info.block_indices.at(local_block_index);

      schur_complement_cholesky_cond.at(global_block_index)
        = cholesky_condition_number(cholesky_block);
    }

  // X_cholesky, Y_cholesky

  const auto num_X_blocks = 2 * num_blocks;
  std::vector<El::BigFloat> X_cholesky_cond(num_X_blocks);
  std::vector<El::BigFloat> Y_cholesky_cond(num_X_blocks);
  ASSERT_EQUAL(X_cholesky.blocks.size(), 2 * block_info.block_indices.size());
  ASSERT_EQUAL(Y_cholesky.blocks.size(), 2 * block_info.block_indices.size());
  for(size_t i = 0; i < block_info.block_indices.size(); ++i)
    {
      const auto block_index = block_info.block_indices.at(i);
      for(const size_t parity : {0, 1})
        {
          const auto X_index = 2 * block_index + parity;
          const auto local_X_index = 2 * i + parity;
          X_cholesky_cond.at(X_index)
            = cholesky_condition_number(X_cholesky.blocks.at(local_X_index));
          Y_cholesky_cond.at(X_index)
            = cholesky_condition_number(Y_cholesky.blocks.at(local_X_index));
        }
    }

  // Synchronize all at rank 0

  El::mpi::Reduce(schur_complement_cholesky_cond.data(), num_blocks,
                  El::mpi::MAX, 0, El::mpi::COMM_WORLD);
  El::mpi::Reduce(X_cholesky_cond.data(), num_X_blocks, El::mpi::MAX, 0,
                  El::mpi::COMM_WORLD);
  El::mpi::Reduce(Y_cholesky_cond.data(), num_X_blocks, El::mpi::MAX, 0,
                  El::mpi::COMM_WORLD);

  if(El::mpi::Rank() == 0)
    {
      for(size_t block_index = 0; block_index < num_blocks; ++block_index)
        {
          if(max_block_cond_number
             < schur_complement_cholesky_cond.at(block_index))
            {
              max_block_cond_number
                = schur_complement_cholesky_cond.at(block_index);
              max_block_cond_number_name = "schur_complement_cholesky.block_"
                                           + std::to_string(block_index);
            }
          for(const size_t parity : {0, 1})
            {
              const auto X_index = 2 * block_index + parity;
              if(max_block_cond_number < X_cholesky_cond.at(X_index))
                {
                  max_block_cond_number = X_cholesky_cond.at(X_index);
                  max_block_cond_number_name = El::BuildString(
                    "X_cholesky.block_", block_index, "_", parity);
                }
              if(max_block_cond_number < Y_cholesky_cond.at(X_index))
                {
                  max_block_cond_number = Y_cholesky_cond.at(X_index);
                  max_block_cond_number_name = El::BuildString(
                    "Y_cholesky.block_", block_index, "_", parity);
                }
            }
        }
    }
}