#pragma once

#include "sdpb_util/assert.hxx"

#include <string>

// How spectrum decides whether a local minimum of det(c-B.y) is a zero.
// See zero_depth.hxx for definitions.
enum class Zero_Test
{
  // Local test on P(x) = det((c-B.y)/prefactor): width of the dip
  // relative to the distance to the nearest other critical point of P.
  curvature,
  // det(c-B.y) at x0 compared with its values at the midpoints between x0
  // and the neighbouring minima: the test of SDPB 3.1.0 and earlier.
  geometric_mean
};

inline Zero_Test to_Zero_Test(const std::string &name)
{
  if(name == "curvature")
    return Zero_Test::curvature;
  if(name == "geometricMean")
    return Zero_Test::geometric_mean;
  RUNTIME_ERROR("Unknown --zeroTest='", name,
                "', allowed values: geometricMean, curvature");
}

inline std::string to_string(const Zero_Test &zero_test)
{
  switch(zero_test)
    {
    case Zero_Test::curvature: return "curvature";
    case Zero_Test::geometric_mean: return "geometricMean";
    }
  RUNTIME_ERROR("Unknown Zero_Test");
}
