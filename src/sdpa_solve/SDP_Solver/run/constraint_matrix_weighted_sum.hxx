#pragma once

#include "sdpa_solve/SDP.hxx"

namespace Sdpb::Sdpa
{
  // result = \sum_p a[p] F_p
  inline void
  constraint_matrix_weighted_sum(const SDP &sdp, const Primal_Dist_Vector &a,
                                 Block_Diagonal_Matrix &result)
  {
    ASSERT_EQUAL(a.Height(), sdp.primal_dimension());
    ASSERT_EQUAL(a.Width(), 1);
    result.set_zero();
    for(size_t i = 0; i < sdp.sdp_blocks_F.size(); i++)
      {
        const auto a_i = a.Get(i, 0);
        const auto &F_i = sdp.sdp_blocks_F.at(i);
        for(size_t block = 0; block < result.blocks.size(); ++block)
          {
            auto x_F_i_block = F_i.blocks.at(block);
            El::Scale(a_i, x_F_i_block);
            result.blocks.at(block) += x_F_i_block;
          }
      }
  }
}
