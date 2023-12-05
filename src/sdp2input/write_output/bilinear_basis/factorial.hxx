#pragma once

#include "sdpb_util/Boost_Float.hxx"

// We implement factorial by hand because boost::math's versions want
// to initialize epsilon in a static constructor.  However, Bigfloat
// has not had its precision set yet, so it ends up dividing by zero
// and crashing.

inline Boost_Float factorial(const int64_t &n)
{
  Boost_Float result(1);
  for(int64_t kk = 2; kk <= n; ++kk)
    {
      result *= kk;
    }
  return result;
}
