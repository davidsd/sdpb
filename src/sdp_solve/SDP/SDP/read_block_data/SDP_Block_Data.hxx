#pragma once

#include "pmp2sdp/Block_File_Format.hxx"
#include "sdp_solve/Block_Info.hxx"

#include <El.hpp>

struct SDP_Block_Data
{
  int block_index_local = -1;
  El::Matrix<El::BigFloat> constraint_matrix{};
  El::Matrix<El::BigFloat> primal_objective_c{};

  std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases{};
  std::array<El::Matrix<El::BigFloat>, 2> bases_blocks{};

  SDP_Block_Data() = default;
  SDP_Block_Data(std::istream &block_stream, Block_File_Format format,
                 size_t block_index_local, const Block_Info &block_info);

  // Allow move and prohibit copy

  SDP_Block_Data(const SDP_Block_Data &other) = delete;
  SDP_Block_Data(SDP_Block_Data &&other) noexcept = default;
  SDP_Block_Data &operator=(const SDP_Block_Data &other) = delete;
  SDP_Block_Data &operator=(SDP_Block_Data &&other) noexcept = default;
};
