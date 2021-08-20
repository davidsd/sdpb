#pragma once

#include <El.hpp>

inline void
fill_weights(const El::Matrix<El::BigFloat> &y, const size_t &max_index,
             const std::vector<El::BigFloat> &normalization,
             std::vector<El::BigFloat> &weights)
{
  // The weight at max_index is determined by the normalization
  // condition dot(norm,weights)=1
  weights.at(max_index) = 1;
  for(size_t block_row(0); block_row != size_t(y.Height()); ++block_row)
    {
      const size_t index(block_row + (block_row < max_index ? 0 : 1));
      weights.at(index) = y(block_row, 0);
      weights.at(max_index) -= weights.at(index) * normalization.at(index);
    }
  weights.at(max_index) /= normalization.at(max_index);
}
