#pragma once

#include <El.hpp>

struct Zeros
{
  std::vector<El::BigFloat> zeros;
  El::Matrix<El::BigFloat> lambda, error;
};
