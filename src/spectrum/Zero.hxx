#pragma once

#include <El.hpp>

struct Zero
{
  El::BigFloat zero;
  El::Matrix<El::BigFloat> lambda, error;
  Zero(const El::BigFloat &z): zero(z) {}
};
