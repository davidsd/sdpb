#pragma once

#include "Zero.hxx"

#include <filesystem>

struct Zeros
{
  std::vector<Zero> zeros;
  El::BigFloat error;
  std::filesystem::path block_path;
};
