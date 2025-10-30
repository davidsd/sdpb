#pragma once

#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.hxx"
#include "sdpb_util/memory_estimates.hxx"

// Q := P^T P
// P: input_height x input_width
// Q: input_width x input_width
struct Bigint_Syrk_Config
{
  El::mpi::Comm shared_memory_comm;
  const size_t group_index;

  const size_t precision;

  const size_t num_nodes;

  // group index (on a node) -> Total height for all blocks of an MPI group
  const std::vector<size_t> &input_height_per_group;

  const size_t input_height;
  const size_t input_width;
  const size_t output_dim;

  const size_t input_split_factor;
  const size_t output_split_factor;

  Bigint_Syrk_Config(const El::mpi::Comm &shared_memory_comm,
                     const size_t group_index, const size_t precision,
                     const size_t num_nodes,
                     const std::vector<size_t> &input_height_per_group,
                     const size_t input_width, const size_t input_split_factor,
                     const size_t output_split_factor)
      : shared_memory_comm(shared_memory_comm),
        group_index(group_index),
        precision(precision),
        num_nodes(num_nodes),
        input_height_per_group(input_height_per_group),
        input_height(std::accumulate(input_height_per_group.begin(),
                                     input_height_per_group.end(), 0)),
        input_width(input_width),
        output_dim(input_width),
        input_split_factor(input_split_factor),
        output_split_factor(output_split_factor)
  {}
  [[nodiscard]] Fmpz_Comb comb() const
  {
    const mp_limb_t bits = precision;
    // Max number of additions during syrk computation
    const size_t max_additions = input_height;
    // Allow negative values
    const int sign = 1;
    return Fmpz_Comb(bits, bits, sign, max_additions);
  }
  [[nodiscard]] size_t num_primes() const { return comb().num_primes; }

  [[nodiscard]] size_t input_window_size() const
  {
    return div_ceil(input_height, input_split_factor) * input_width
           * num_primes();
  }
  [[nodiscard]] size_t num_input_windows() const
  {
    return output_split_factor == 1 ? 1 : 2;
  }
  [[nodiscard]] size_t output_window_size() const
  {
    const size_t dim = div_ceil(output_dim, output_split_factor);
    return dim * dim * num_primes();
  }
  // Copied from BigInt_Shared_Memory_Syrk_Context.cxx
  // Estimate total MPI buffer sizes on a node for restore_and_reduce()
  [[nodiscard]] size_t get_reduce_scatter_buffer_size() const
  {
    // no reduce-scatter => no RAM needed
    if(num_nodes == 1)
      return 0;

    const size_t dim = div_ceil(output_dim, output_split_factor);

    // Total number of elements in the output submatrix
    // If split factor == 1, we need only half of the matrix.
    // Otherwise, we need the whole matrix for off-diagonal blocks
    const size_t num_elements
      = output_split_factor == 1 ? dim * (dim + 1) / 2 : dim * dim;

    // Each rank needs ~(num_elements/num_ranks) buffer for MPI_Send
    // and the same for MPI_Recv
    return 2 * div_ceil(num_elements, num_nodes);
  }

  // Total memory estimates

  // Local memory allocated on a node
  [[nodiscard]] size_t node_local_bytes() const
  {
    return bigfloat_bytes(precision) * get_reduce_scatter_buffer_size();
  }
  // Local memory allocated on a node
  [[nodiscard]] size_t node_shmem_bytes() const
  {
    // Input and output residue windows
    return sizeof(double)
           * (num_input_windows() * input_window_size()
              + output_window_size());
  }
  // All memory allocated on a node
  [[nodiscard]] size_t node_total_bytes() const
  {
    return node_local_bytes() + node_shmem_bytes();
  }

  // Helper functions. TODO move elsewhere and reuse?
private:
  static size_t div_ceil(const size_t a, const size_t b)
  {
    return a / b + (a % b == 0 ? 0 : 1);
  }
};