#pragma once

#include "Environment.hxx"
#include "Verbosity.hxx"

size_t
get_max_shared_memory_bytes(size_t nonshared_memory_required_per_node_bytes,
                            const Environment &env, Verbosity verbosity);

size_t bigfloat_bytes();

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
