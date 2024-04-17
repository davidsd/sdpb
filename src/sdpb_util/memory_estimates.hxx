#pragma once

#include "Environment.hxx"
#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdp_solve/SDP_Solver.hxx"

size_t
get_max_shared_memory_bytes(size_t nonshared_memory_required_per_node_bytes,
                            const Environment &env, Verbosity verbosity);

size_t get_matrix_size_local(const Block_Diagonal_Matrix &X);
size_t get_A_X_size_local(const Block_Info &block_info, const SDP &sdp);
size_t get_schur_complement_size_local(const Block_Info &block_info);
size_t get_B_size_local(const SDP &sdp);
size_t get_Q_size_local(const SDP &sdp);
size_t get_SDP_size_local(const SDP &sdp);

size_t bigfloat_bytes();

size_t get_heap_allocated_bytes(const El::BigFloat &f);
size_t get_heap_allocated_bytes(const Block_Diagonal_Matrix &m);
size_t get_heap_allocated_bytes(const Block_Matrix &m);
size_t get_heap_allocated_bytes(const Block_Vector &v);
size_t get_heap_allocated_bytes(const El::DistMatrix<El::BigFloat> &m);
size_t get_heap_allocated_bytes(const El::Matrix<El::BigFloat> &m);
size_t get_heap_allocated_bytes(const SDP &sdp);
size_t get_heap_allocated_bytes(const SDP_Solver &solver);
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
