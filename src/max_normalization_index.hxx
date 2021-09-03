#pragma once

#include <El.hpp>

inline size_t max_normalization_index(const std::vector<El::BigFloat> &normalization)
{
  auto max_normalization(normalization.begin());
  for(auto n(normalization.begin()); n != normalization.end(); ++n)
    {
      if(Abs(*n) > Abs(*max_normalization))
        {
          max_normalization = n;
        }
    }
  return std::distance(normalization.begin(), max_normalization);
}
