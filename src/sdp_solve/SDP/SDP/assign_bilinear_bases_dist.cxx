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
      for(int64_t row = 0; row < dist.LocalHeight(); ++row)
        {
          El::Int global_row(dist.GlobalRow(row));
          for(int64_t column = 0; column < dist.LocalWidth(); ++column)
            {
              El::Int global_column(dist.GlobalCol(column));
              dist.SetLocal(row, column, local(global_row, global_column));
            }
        }
    }
}
