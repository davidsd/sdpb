#pragma once

#include "sdp_solve/Block_Info.hxx"

void set_bilinear_bases_block_local(
  const El::Matrix<El::BigFloat> &bilinear_base_local,
  El::Matrix<El::BigFloat> &bases_block_local);

void set_bilinear_bases_block(
  const El::Matrix<El::BigFloat> &bilinear_base_local,
  El::DistMatrix<El::BigFloat> &bases_block);

void set_bases_blocks(
  const Block_Info &block_info,
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  const El::Grid &grid);
