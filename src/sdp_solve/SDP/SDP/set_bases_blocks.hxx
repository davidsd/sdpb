#pragma once

#include "sdp_solve/Block_Info.hxx"

void set_bases_blocks(
  const Block_Info &block_info,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  const El::Grid &grid);
