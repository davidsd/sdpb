#pragma once

#include "sdpa_solve/Block_Info.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.hxx"

namespace Sdpb::Sdpa
{
  struct Initialize_P_Config
  {
    El::mpi::Comm shared_memory_comm;
    const size_t precision;
    const Block_Mapping block_mapping;
    const size_t node_index;
    const size_t group_index;
    const std::vector<size_t> node_block_dims;
    const size_t max_block_dim;
    const std::vector<size_t> max_group_block_dims;
    const size_t sum_dim_squared;
    const size_t primal_dimension;
    const size_t split_factor;

    Initialize_P_Config(const Block_Info &block_info, size_t split_factor);
    Initialize_P_Config(El::mpi::Comm shared_memory_comm, size_t precision,
                        const Block_Mapping &block_mapping, size_t node_index,
                        size_t group_index,
                        const std::vector<size_t> &block_dimensions,
                        size_t primal_dimension, size_t split_factor);

    const std::vector<Block_Location> &local_block_locations() const;
    const std::vector<Block_Location> &node_block_locations() const;

    [[nodiscard]] Fmpz_Comb comb() const;
    [[nodiscard]] size_t num_primes() const;
    // Total number of elements in triangular block-diagonal matrix L
    [[nodiscard]] size_t L_size() const;
    // Total number of matrix elements in a buffer for processing
    // vector<Block_Diagonal_Matrix> G_p (transformations of F_p).
    // We process primal_dimension_step elements of this vector at once.
    [[nodiscard]] size_t G_size() const;
    [[nodiscard]] size_t P_size() const;
    [[nodiscard]] size_t L_normalizer_size() const;
    [[nodiscard]] size_t G_normalizer_size() const;
    [[nodiscard]] size_t L_window_size() const;
    [[nodiscard]] size_t G_window_size() const;

    // Total memory estimates on a node

    [[nodiscard]] size_t L_X_inv_window_bytes() const;
    [[nodiscard]] size_t L_Y_window_bytes() const;
    [[nodiscard]] size_t G_window_bytes() const;

    [[nodiscard]] size_t node_P_bytes() const;
    [[nodiscard]] size_t node_L_bytes() const;
    [[nodiscard]] size_t node_L_normalizer_bytes() const;
    [[nodiscard]] size_t node_G_bytes() const;
    [[nodiscard]] size_t node_G_Normalizer_bytes() const;
    // Temporary allocations inside normalize_and_shift()
    [[nodiscard]] size_t node_normalize_and_shift_L_bytes() const;
    [[nodiscard]] size_t node_normalize_and_shift_G_bytes() const;
    [[nodiscard]] size_t node_reshape_bytes() const;

    [[nodiscard]] size_t max_primal_dimension_step() const;
  };
}