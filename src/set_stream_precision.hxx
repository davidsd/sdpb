#pragma once

#include <El.hpp>
#include <cmath>

inline void set_stream_precision(std::ofstream &ofs)
{
  ofs.precision(
    std::ceil(El::gmp::Precision() * std::log(2.0) / std::log(10.0)));
}
