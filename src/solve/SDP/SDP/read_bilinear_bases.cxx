#include "read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>

void read_bilinear_bases(
  const boost::filesystem::path &sdp_directory, const El::Grid &grid,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases_dist)
{
  boost::filesystem::ifstream bilinear_stream(sdp_directory
                                              / "bilinear_bases");
  if(!bilinear_stream.good())
    {
      throw std::runtime_error("Could not open '"
                               + (sdp_directory / "bilinear_bases").string()
                               + "'");
    }

  size_t size;
  bilinear_stream >> size;
  bilinear_bases_local.reserve(size);
  bilinear_bases_dist.reserve(size);
  for(size_t m = 0; m < size; ++m)
    {
      size_t height, width;
      bilinear_stream >> height >> width;

      bilinear_bases_local.emplace_back(height, width);
      auto &local(bilinear_bases_local.back());
      for(size_t row = 0; row < height; ++row)
        for(size_t column = 0; column < width; ++column)
          {
            bilinear_stream >> local(row, column);
          }

      bilinear_bases_dist.emplace_back(height, width, grid);
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
