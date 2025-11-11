#pragma once

#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.hxx"
#include "sdpb_util/memory_estimates.hxx"

// Q := P^T P
// P: input_height x input_width
// Q: input_width x input_width
struct Bigint_Syrk_Config
{
  El::mpi::Comm shared_memory_comm;
  const size_t precision;
  const size_t num_nodes;

  // group index (on a node) -> Total height for all blocks of an MPI group
  const std::vector<size_t> input_height_per_group;

  const size_t input_height;
  const size_t input_width;
  const size_t output_dim;

  const size_t input_split_factor;
  const size_t output_split_factor;

  const std::vector<size_t> input_window_height_per_group_per_prime;
  // const size_t input_window_width;

  Bigint_Syrk_Config(const El::mpi::Comm &shared_memory_comm,
                     const size_t precision, const size_t num_nodes,
                     const std::vector<size_t> &input_height_per_group,
                     const size_t input_width, const size_t input_split_factor,
                     const size_t output_split_factor);
  [[nodiscard]] Fmpz_Comb comb() const;
  [[nodiscard]] size_t num_primes() const;

  [[nodiscard]] size_t input_window_height() const;
  // Input and output windows have the same width
  [[nodiscard]] size_t window_width() const;
  [[nodiscard]] size_t input_window_size() const;
  [[nodiscard]] size_t num_input_windows() const;
  [[nodiscard]] size_t output_window_size() const;
  // Copied from BigInt_Shared_Memory_Syrk_Context.cxx
  // Estimate total MPI buffer sizes on a node for restore_and_reduce()
  [[nodiscard]] size_t get_reduce_scatter_buffer_size() const;

  // Total memory estimates

  // Local memory allocated on a node
  [[nodiscard]] size_t node_local_bytes() const;
  // Local memory allocated on a node
  [[nodiscard]] size_t node_shmem_bytes() const;
  // All memory allocated on a node
  [[nodiscard]] size_t node_total_bytes() const;
};