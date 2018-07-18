#include "../../Block_Matrix.hxx"

#include <boost/filesystem.hpp>

void read_free_var_matrix(const boost::filesystem::path &sdp_directory,
                          const std::vector<size_t> &block_indices,
                          const El::Grid &grid, Block_Matrix &free_var_matrix)
{
  for(auto &block_index : block_indices)
    {
      boost::filesystem::ifstream free_var_matrix_stream(
        sdp_directory / ("free_var_matrix." + std::to_string(block_index)));

      size_t height, width;
      free_var_matrix_stream >> height >> width;

      El::Matrix<El::BigFloat> temp(height, width);
      for(size_t row = 0; row < height; ++row)
        for(size_t column = 0; column < width; ++column)
          {
            free_var_matrix_stream >> temp(row, column);
          }

      free_var_matrix.blocks.emplace_back(height, width, grid);
      auto &block(free_var_matrix.blocks.back());
      size_t local_height(block.LocalHeight()),
        local_width(block.LocalWidth());
      for(size_t row = 0; row < local_height; ++row)
        {
          size_t global_row(block.GlobalRow(row));
          for(size_t column = 0; column < local_width; ++column)
            {
              size_t global_column(block.GlobalCol(column));
              block.SetLocal(row, column, temp(global_row, global_column));
            }
        }
    }
}
