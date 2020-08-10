#include "../../SDP.hxx"
#include "../../read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_bilinear_bases(
  const boost::filesystem::path &sdp_directory, const Block_Info &block_info,
  const El::Grid &grid,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local,
  std::vector<El::DistMatrix<El::BigFloat>> &bilinear_bases_dist)
{
  auto &block_indices(block_info.block_indices);
  bilinear_bases_local.resize(2 * block_indices.size());

  for(size_t file_rank(0); file_rank < block_info.file_num_procs; ++file_rank)
    {
      const boost::filesystem::path bilinear_path(
        sdp_directory / ("bilinear_bases." + std::to_string(file_rank)));
      boost::filesystem::ifstream bilinear_stream(bilinear_path);
      if(!bilinear_stream.good())
        {
          throw std::runtime_error("Could not open '" + bilinear_path.string()
                                   + "'");
        }
      size_t file_num_bases;
      bilinear_stream >> file_num_bases;
      if(!bilinear_stream.good())
        {
          throw std::runtime_error("Corrupted or empty file: "
                                   + bilinear_path.string());
        }

      // The index of bilinear_bases in each file is described by
      // file_block_indices.  However, block_indices is not in
      // numerical order.  So we have to take care when placing blocks
      for(size_t block = 0; block < file_num_bases; ++block)
        for(size_t parity = 0; parity < 2; ++parity)
          {
            size_t height, width;
            bilinear_stream >> height >> width;
            if(!bilinear_stream.good())
              {
                throw std::runtime_error("Corrupted header in file: "
                                         + bilinear_path.string());
              }

            auto block_iter(std::find(
              block_indices.begin(), block_indices.end(),
              block_info.file_block_indices.at(file_rank).at(block)));

            if(block_iter != block_indices.end())
              {
                size_t block_index(
                  2 * std::distance(block_indices.begin(), block_iter)
                  + parity);
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
      if(!bilinear_stream.good())
        {
          throw std::runtime_error("Corrupted data in file: "
                                   + bilinear_path.string());
        }
    }

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
