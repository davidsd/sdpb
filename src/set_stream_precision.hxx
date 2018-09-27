#pragma once

#include <El.hpp>
#include <cmath>
#include <iostream>

inline void set_stream_precision(std::ostream &os)
{
  // This is the same formula as used by std::numeric_limits::max_digits10
  os.precision(std::ceil(El::gmp::Precision() * std::log10(2.0)) + 1);
}
