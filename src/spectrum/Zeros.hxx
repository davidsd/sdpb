#pragma once

#include <El.hpp>

#include <deque>

struct Zeros
{
  std::deque<El::BigFloat> zeros;
  El::Matrix<El::BigFloat> lambda, error;
};
