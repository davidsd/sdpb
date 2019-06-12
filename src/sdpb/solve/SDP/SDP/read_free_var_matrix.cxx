#include "../../Block_Matrix.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_free_var_matrix(const boost::filesystem::path &sdp_directory,
                          const std::vector<size_t> &block_indices,
                          const El::Grid &grid, Block_Matrix &B)
{
  B.blocks.reserve(block_indices.size());
  for(auto &block_index : block_indices)
    {
      boost::filesystem::ifstream B_stream(
        sdp_directory / ("free_var_matrix." + std::to_string(block_index)));
      if(!B_stream.good())
        {
          throw std::runtime_error(
            "Could not open '"
            + (sdp_directory
               / ("free_var_matrix." + std::to_string(block_index)))
                .string()
            + "'");
        }

      size_t height, width;
      B_stream >> height >> width;
      B.blocks.emplace_back(height, width, grid);
      auto &block(B.blocks.back());
      for(size_t row = 0; row < height; ++row)
        for(size_t column = 0; column < width; ++column)
          {
            El::BigFloat input_num;
            B_stream >> input_num;
            if(block.IsLocal(row, column))
              {
                block.SetLocal(block.LocalRow(row), block.LocalCol(column),
                               input_num);
              }
          }
    }
}
