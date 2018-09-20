#include "../../../read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>

void read_bilinear_bases(
  const boost::filesystem::path &sdp_directory,
  const std::vector<size_t> &block_indices, const El::Grid &grid,
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
  bilinear_bases_local.resize(2 * block_indices.size());
  // bilinear_bases is in numerical order, but block_indices is not.
  // So we have to take care when placing blocks
  for(size_t block = 0; block < size; ++block)
    for(size_t parity = 0; parity < 2; ++parity)
      {
        size_t height, width;
        bilinear_stream >> height >> width;
        auto block_iter(
          std::find(block_indices.begin(), block_indices.end(), block));
        if(block_iter != block_indices.end())
          {
            size_t block_index(
              2 * std::distance(block_indices.begin(), block_iter) + parity);
            auto &local(bilinear_bases_local.at(block_index));
            local.Resize(height, width);
            for(size_t row = 0; row < height; ++row)
              for(size_t column = 0; column < width; ++column)
                {
                  bilinear_stream >> local(row, column);
                }
          }
        else
          {
            // Add one to get the initial newline after 'width'.
            for(size_t line = 0; line < height * width + 1; ++line)
              {
                bilinear_stream.ignore(
                  std::numeric_limits<std::streamsize>::max(), '\n');
              }
          }
      }

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
