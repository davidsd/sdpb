#pragma once

#include <El.hpp>

struct Zero
{
  El::BigFloat zero;
  El::Matrix<El::BigFloat> lambda;
  Zero(const El::BigFloat &z) : zero(z) {}
};
