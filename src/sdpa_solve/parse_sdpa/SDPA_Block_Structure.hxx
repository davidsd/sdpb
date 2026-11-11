#pragma once

#include "sdpb_util/assert.hxx"

#include <vector>

namespace Sdpb::Sdpa
{
  struct SDPA_Block_Structure
  {
    const size_t m_dim;
    // Can contain negative numbers which denote diagonal blocks.
    const std::vector<int> input_block_sizes;
    // All sizes are positive. Each input diagonal block is split into many 1x1 solver blocks.
    // TODO: shall we use El::SparseMatrix instead?
    std::vector<size_t> solver_block_sizes;
    size_t num_solver_blocks;

  private:
    std::vector<size_t> num_prev_blocks;

  public:
    SDPA_Block_Structure(const size_t m, const size_t nblocks,
                         const std::vector<int> &block_sizes)
        : m_dim(m), input_block_sizes(block_sizes)
    {
      ASSERT_EQUAL(nblocks, block_sizes.size());

      size_t prev_blocks = 0;
      for(const auto &block_size : block_sizes)
        {
          num_prev_blocks.push_back(prev_blocks);
          ASSERT(block_size != 0);
          if(block_size > 0)
            {
              prev_blocks += 1;
              solver_block_sizes.push_back(block_size);
            }
          // if block_size is negative, this means a diagonal block.
          // We split it into many 1x1 blocks
          else
            {
              prev_blocks += -block_size;
              solver_block_sizes.insert(solver_block_sizes.end(), -block_size,
                                        1);
            }
        }
      num_solver_blocks = prev_blocks;
    }

    void get_element_position(const size_t input_block_index,
                              const size_t input_i, const size_t input_j,
                              // output parameters:
                              size_t &solver_block_index, size_t &i,
                              size_t &j) const
    {
      const auto &block_size = input_block_sizes.at(input_block_index - 1);
      ASSERT(block_size != 0);
      ASSERT(input_i > 0);
      ASSERT(input_i > 0);
      const auto prev_blocks = num_prev_blocks.at(input_block_index - 1);
      if(block_size > 0)
        {
          solver_block_index = prev_blocks;
          i = input_i - 1;
          j = input_j - 1;
        }
      else
        {
          // We represent a diagonal SDP block as a collection of 1x1 blocks.
          ASSERT_EQUAL(input_i, input_j,
                       "Non-diagonal element found in diagonal SDP block!",
                       DEBUG_STRING(input_block_index),
                       DEBUG_STRING(block_size));
          solver_block_index = prev_blocks + input_i - 1;
          i = 0;
          j = 0;
        }
    }
  };
}
