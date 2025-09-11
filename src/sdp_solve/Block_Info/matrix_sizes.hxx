#pragma once

#include <El.hpp>
#include <vector>

inline size_t
get_schur_block_size(const std::vector<size_t> &dimensions,
                     const std::vector<size_t> &num_points, const size_t index)
{
  return num_points.at(index) * dimensions.at(index)
         * (dimensions.at(index) + 1) / 2;
}
inline std::vector<size_t>
schur_block_sizes(const std::vector<size_t> &dimensions,
                  const std::vector<size_t> &num_points)
{
  std::vector<size_t> result(num_points.size());
  for(size_t index(0); index != num_points.size(); ++index)
    {
      result[index] = get_schur_block_size(dimensions, num_points, index);
    }
  return result;
}
inline size_t
get_bilinear_pairing_block_size(const std::vector<size_t> &dimensions,
                                const std::vector<size_t> &num_points,
                                const size_t index, const size_t parity)
{
  ASSERT(parity == 0 || parity == 1, DEBUG_STRING(parity));
  return num_points.at(index) * dimensions.at(index);
}
inline std::vector<size_t>
bilinear_pairing_block_sizes(const std::vector<size_t> &dimensions,
                             const std::vector<size_t> &num_points)
{
  std::vector<size_t> result(2 * num_points.size());
  for(size_t index(0); index != num_points.size(); ++index)
    {
      result[2 * index]
        = get_bilinear_pairing_block_size(dimensions, num_points, index, 0);
      result[2 * index + 1]
        = get_bilinear_pairing_block_size(dimensions, num_points, index, 1);
    }
  return result;
}
inline size_t
get_psd_matrix_block_size(const std::vector<size_t> &dimensions,
                          const std::vector<size_t> &num_points,
                          const size_t index, const size_t parity)
{
  // Need to round down (num_points+1)/2 before multiplying by
  // dim, since dim could be 2.
  const size_t even = dimensions.at(index) * ((num_points.at(index) + 1) / 2);
  if(parity == 0)
    return even;
  if(parity == 1)
    return dimensions.at(index) * num_points.at(index) - even;
  RUNTIME_ERROR("parity should be 0 or 1", DEBUG_STRING(parity));
}
inline std::vector<size_t>
psd_matrix_block_sizes(const std::vector<size_t> &dimensions,
                       const std::vector<size_t> &num_points)
{
  std::vector<size_t> result(2 * num_points.size());
  for(size_t index(0); index != num_points.size(); ++index)
    {
      for(const size_t parity : {0, 1})
        {
          result[2 * index + parity]
            = get_psd_matrix_block_size(dimensions, num_points, index, parity);
        }
    }
  return result;
}
