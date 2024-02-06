#pragma once

#include "outer_limits/Function.hxx"

#include <El.hpp>

#include <vector>

struct Functions_File_Parse_Result
{
  std::vector<El::BigFloat> objective;
  std::vector<El::BigFloat> normalization;
  std::vector<std::vector<std::vector<std::vector<Function>>>> functions;

  Functions_File_Parse_Result() = default;

  // Allow moving and prevent accidential copying

  Functions_File_Parse_Result(const Functions_File_Parse_Result &other)
    = delete;
  Functions_File_Parse_Result(Functions_File_Parse_Result &&other) noexcept
    = default;
  Functions_File_Parse_Result &
  operator=(const Functions_File_Parse_Result &other)
    = delete;
  Functions_File_Parse_Result &
  operator=(Functions_File_Parse_Result &&other) noexcept
    = default;
};
