#pragma once

#include "Block_Info.hxx"
#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>

namespace Sdpb::Sdpa
{
  // c_1..c_m, x_1..x_m etc.
  using Primal_Dist_Vector = El::DistMatrix<El::BigFloat, El::STAR, El::STAR>;

  // minimize sum_{i=1}^m c_i x_i
  // s.t. X = sum_{i=1}^m F_i x_i - F_0 is positive semidefinite.
  // Matrices F_i have block diagonal structure
  struct SDP
  {
    // matrices F_1, ... F_m
    std::vector<Block_Diagonal_Matrix> sdp_blocks_F;

    Block_Diagonal_Matrix sdp_block_F_0;

    // c_1..c_m, a vector of length m
    // It is duplicated amongst all the blocks.
    Primal_Dist_Vector primal_objective_c;

    SDP(const std::filesystem::path &sdp_path, const Block_Info &block_info,
        const El::Grid &grid, Timers &timers);

    size_t primal_dimension() const;

  private:
    void validate(const Block_Info &block_info) const noexcept(false);
  };
}
