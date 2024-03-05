#pragma once

#include "sdpb_util/Boost_Float.hxx"

#include <numeric>

template <typename T>
Boost_Float accumulate_over_others(
  const std::vector<Boost_Float> &sorted_poles,
  const std::pair<std::vector<Boost_Float>::const_iterator,
                  std::vector<Boost_Float>::const_iterator> &equal_range,
  const Boost_Float &start_value, T op)
{
  Boost_Float first_part(
    std::accumulate(sorted_poles.begin(), equal_range.first, start_value, op));
  return std::accumulate(equal_range.second, sorted_poles.end(), first_part,
                         op);
}
