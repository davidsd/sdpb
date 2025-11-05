#pragma once

#include "Environment.hxx"
#include "Verbosity.hxx"
#include "assert.hxx"

#include <algorithm>

size_t
get_max_shared_memory_bytes(size_t nonshared_memory_required_per_node_bytes,
                            const Environment &env, Verbosity verbosity);

size_t bigfloat_bytes();
size_t bigfloat_bytes(mp_bitcnt_t precision);

size_t get_heap_allocated_bytes(const El::BigFloat &f);
size_t get_heap_allocated_bytes(const El::AbstractDistMatrix<El::BigFloat> &m);
size_t get_heap_allocated_bytes(const El::Matrix<El::BigFloat> &m);
template <class T> size_t get_heap_allocated_bytes(const std::vector<T> &vec)
{
  size_t res = 0;
  for(const auto &element : vec)
    res += sizeof(element) + get_heap_allocated_bytes(element);
  res += sizeof(T) * (vec.capacity() - vec.size());
  return res;
}
template <class T, std::size_t N>
size_t get_heap_allocated_bytes(const std::array<T, N> &arr)
{
  size_t res = 0;
  for(const auto &element : arr)
    res += get_heap_allocated_bytes(element);
  return res;
}
template <class T> size_t get_allocated_bytes(const T &value)
{
  return sizeof(value) + get_heap_allocated_bytes(value);
}

// Sum bytes for all ranks on a node and print from rank=0
void print_allocation_message_per_node(const Environment &env,
                                       const std::string &name, size_t bytes);

// Extra memory required to compute L^{-1} X, via
// El::Trsm(LEFT, LOWER, NORMAL, NON_UNIT, 1, L, X);
// L ~ height x height, lower triangular
// X ~ height x width
// See https://gitlab.com/bootstrapcollaboration/elemental/-/blob/3ee9f404e632d9d0a4751e9b6717da51e59b5697/src/blas_like/level3/Trsm.cpp#L119-L122
//
// NB: if grid_size (=MPI group size) is a prime number (e.g. 11, 13, or 17),
// then the grid has 1D shape (e.g. 11x1, 13x1, or 17x1).
// In this case the terms proportional to grid_height will lead to
// very high memory consumption compared to more balanced grids (e.g. 4x3 or 4x4).
// To avoid this issue, set e.g. --procGranularity 2 or --procGranularity 4.
size_t get_trsm_bytes(int height, int width, int grid_height, int grid_width);

// Extra memory required to compute X L, via
// El::Trmm(RIGHT, LOWER, NORMAL, NON_UNIT, 1, L, X)
// X ~ height x width
// L ~ width x width, lower triangular
// See https://gitlab.com/bootstrapcollaboration/elemental/-/blob/3ee9f404e632d9d0a4751e9b6717da51e59b5697/src/blas_like/level3/Trmm.cpp#L78
//
// NB: if grid_size (=MPI group size) is a prime number (e.g. 11, 13, or 17),
// then the grid has 1D shape (e.g. 11x1, 13x1, or 17x1).
// In this case the terms proportional to grid_height will lead to
// very high memory consumption compared to more balanced grids (e.g. 4x3 or 4x4).
// To avoid this issue, set e.g. --procGranularity 2 or --procGranularity 4.
size_t get_trmm_bytes(int height, int width, int grid_height, int grid_width);
