#pragma once

#include <El.hpp>

#include <vector>
#include <map>

struct Function_Block
{
  std::vector<El::BigFloat> points;
  std::vector<std::map<El::BigFloat,El::BigFloat>> functions;
};
