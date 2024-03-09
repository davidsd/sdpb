#pragma once

#include "assert.hxx"

#include <El.hpp>

#include <vector>

// Split as even as possible,
// e.g. for split_factor=3 the range [0,10)
// is split into [0,4), [4,7), [7,10)
std::vector<El::Range<El::Int>> inline split_range(
  const El::Range<El::Int> &range, size_t split_factor)
{
  ASSERT(range.beg >= 0 && range.end > range.beg, DEBUG_STRING(range.beg),
         DEBUG_STRING(range.end));

  const size_t total_size = range.end - range.beg;
  // We cannot split into more than total_size ranges
  split_factor = std::min(split_factor, total_size);

  std::vector<El::Range<El::Int>> ranges(split_factor);
  El::Int start = range.beg;
  for(size_t n = 0; n < split_factor; ++n)
    {
      El::Int dim = total_size / split_factor;
      if(n < total_size % split_factor)
        dim++;
      ASSERT(dim > 0, DEBUG_STRING(total_size), DEBUG_STRING(split_factor),
             DEBUG_STRING(n));

      ranges.at(n) = {start, start + dim};
      start += dim;
    }
  return ranges;
}