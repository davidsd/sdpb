#pragma once

#include <El.hpp>

#include <filesystem>

inline void
set_element(El::DistMatrix<El::BigFloat> &block, const int64_t &row,
            const int64_t &column, const El::BigFloat &element)
{
  if(block.IsLocal(row, column))
    {
      block.SetLocal(block.LocalRow(row), block.LocalCol(column),
                     El::BigFloat(element));
    }
}

inline void set_element(El::Matrix<El::BigFloat> &block, const int64_t &row,
                        const int64_t &column, const El::BigFloat &element)
{
  block(row, column) = element;
}

template <typename Matrix>
void read_text_block(Matrix &block, const std::filesystem::path &block_path)
{
  std::ifstream block_stream(block_path);
  if(!block_stream)
    {
      RUNTIME_ERROR("Unable to open checkpoint file: ", block_path);
    }
  int64_t file_height, file_width;
  block_stream >> file_height >> file_width;
  if(!block_stream.good())
    {
      RUNTIME_ERROR("Corrupted header in file: ", block_path);
    }
  ASSERT(file_height == block.Height() && file_width == block.Width(),
         "Incompatible checkpoint file: ", block_path,
         ":  Expected dimensions (", block.Height(), ",", block.Width(),
         "), but found (", file_height, ",", file_width, ")");

  std::string element;
  for(int64_t row = 0; row < file_height; ++row)
    for(int64_t column = 0; column < file_width; ++column)
      {
        block_stream >> element;
        set_element(block, row, column, El::BigFloat(element));
      }
  ASSERT(block_stream.good(), "Corrupted data in file: ", block_path);
}

template <typename Matrix>
void read_text_block(Matrix &block,
                     const std::filesystem::path &block_directory,
                     const std::string &prefix, const size_t &block_index)
{
  read_text_block(block, block_directory
                           / (prefix + std::to_string(block_index) + ".txt"));
}
