#pragma once

#include <El.hpp>

#include <vector>

struct Block_Data_Parse_Result
{
  El::Matrix<El::BigFloat> B;
  std::vector<El::BigFloat> c;
  El::Matrix<El::BigFloat> bilinear_bases_even;
  El::Matrix<El::BigFloat> bilinear_bases_odd;

  Block_Data_Parse_Result() = default;

  void clear()
  {
    c.clear();
    B.Empty();
    bilinear_bases_even.Empty();
    bilinear_bases_odd.Empty();
  }

  // Allow move and prohibit copy

  Block_Data_Parse_Result(const Block_Data_Parse_Result &other) = delete;
  Block_Data_Parse_Result(Block_Data_Parse_Result &&other) noexcept = default;
  Block_Data_Parse_Result &operator=(const Block_Data_Parse_Result &other)
    = delete;
  Block_Data_Parse_Result &operator=(Block_Data_Parse_Result &&other) noexcept
    = default;
};
