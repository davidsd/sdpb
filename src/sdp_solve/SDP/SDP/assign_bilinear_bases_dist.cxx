#include "sdpb_util/copy_matrix.hxx"

#include <El.hpp>

void assign_bilinear_bases_dist(
  const std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  const El::Grid &grid,
  std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases_dist)
{
  bilinear_bases_dist.reserve(bilinear_bases_local.size());
  for(auto &local : bilinear_bases_local)
    {
      bilinear_bases_dist.emplace_back(local.Height(), local.Width(), grid);
      auto &dist(bilinear_bases_dist.back());
      copy_matrix(local, dist);
    }
}
