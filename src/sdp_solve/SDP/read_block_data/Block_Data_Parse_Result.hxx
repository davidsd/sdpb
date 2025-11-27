#pragma once

#include <El.hpp>
#include <optional>

#include <vector>

struct Block_Data_Parse_Result
{
  El::Matrix<El::BigFloat> B;
  std::vector<El::BigFloat> c;
  std::optional<std::vector<El::BigFloat>> preconditioning_values;
  El::Matrix<El::BigFloat> bilinear_bases_even;
  El::Matrix<El::BigFloat> bilinear_bases_odd;

  Block_Data_Parse_Result() = default;

  void clear()
  {
    c.clear();
    preconditioning_values.reset();
    B.Empty();
    bilinear_bases_even.Empty();
    bilinear_bases_odd.Empty();
  }

  // Allow move and prohibit copy

  Block_Data_Parse_Result(const Block_Data_Parse_Result &other) = delete;
  Block_Data_Parse_Result(Block_Data_Parse_Result &&other) = default;
  Block_Data_Parse_Result &operator=(const Block_Data_Parse_Result &other)
    = delete;
  Block_Data_Parse_Result &operator=(Block_Data_Parse_Result &&other)
    = default;
};
