#pragma once

#include <El.hpp>
#include <vector>

struct Function
{
  El::BigFloat max_delta, infinity_value;
  std::vector<El::BigFloat> chebyshev_coeffs;

  El::BigFloat eval(const El::BigFloat &infinity,
                    const El::BigFloat &x) const;
};
