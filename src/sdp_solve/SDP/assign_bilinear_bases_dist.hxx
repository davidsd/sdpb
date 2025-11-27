#pragma once

#include <El.hpp>

void assign_bilinear_bases_dist(
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  const El::Grid &grid,
  std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases_dist);
