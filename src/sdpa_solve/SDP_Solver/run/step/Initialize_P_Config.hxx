#pragma once

#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.hxx"
#include "sdpb_util/memory_estimates.hxx"
#include "sdpb_util/block_mapping/Block_Mapping.hxx"

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
    const size_t sum_dim_squared;
    const size_t primal_dimension;
    // TODO find better name
    // TODO: use split_factor instead?
    const size_t primal_dimension_step;

    Initialize_P_Config(const El::mpi::Comm shared_memory_comm,
                        const size_t precision,
                        const Block_Mapping &block_mapping,
                        const size_t node_index, const size_t group_index,
                        const std::vector<size_t> &block_dimensions,
                        const size_t primal_dimension,
                        const size_t primal_dimension_step)
        : shared_memory_comm(shared_memory_comm),
          precision(precision),
          block_mapping(block_mapping),
          node_index(node_index),
          group_index(group_index),
          node_block_dims(
            get_node_block_dims(node_block_locations(), block_dimensions)),
          max_block_dim(get_max(node_block_dims)),
          sum_dim_squared(sum_squares(node_block_dims)),
          primal_dimension(primal_dimension),
          primal_dimension_step(primal_dimension_step)
    {}

    const std::vector<Block_Location> &local_block_locations() const
    {
      return block_mapping.node_group_block_locations.at(node_index)
        .at(group_index);
    }
    const std::vector<Block_Location> &node_block_locations() const
    {
      return block_mapping.node_block_locations.at(node_index);
    }

    [[nodiscard]] Fmpz_Comb comb() const
    {
      const mp_limb_t bits = precision;
      // Max number of additions during trmm computation.
      const size_t max_additions = max_block_dim;
      // Allow negative values
      constexpr int sign = 1;
      return Fmpz_Comb(bits, bits, sign, max_additions);
    }
    [[nodiscard]] size_t num_primes() const { return comb().num_primes; }
    // Total number of elements in triangular block-diagonal matrix L
    [[nodiscard]] size_t L_size() const { return sum_dim_squared; }
    // Total number of matrix elements in a buffer for processing
    // vector<Block_Diagonal_Matrix> F_p.
    // We process primal_dimension_step elements of this vector at once.
    [[nodiscard]] size_t F_size() const
    {
      return sum_dim_squared * primal_dimension_step;
    }
    [[nodiscard]] size_t P_size() const
    {
      return sum_dim_squared * primal_dimension;
    }
    [[nodiscard]] size_t L_window_size() const
    {
      return L_size() * num_primes();
    }
    [[nodiscard]] size_t F_window_size() const
    {
      return F_size() * num_primes();
    }

    // Total memory estimates

    // Local memory allocated on a node
    [[nodiscard]] size_t node_local_bytes() const
    {
      // L_X_inv, L_Y, F_{p_min..p_max} (also L_X_inv F_p etc.), P
      return bigfloat_bytes(precision) * (2 * L_size() + F_size() + P_size());
    }
    // Local memory allocated on a node
    [[nodiscard]] size_t node_shmem_bytes() const
    {
      // residues for L_X_inv, L_Y, F_p
      return sizeof(double) * (2 * L_window_size() + F_window_size());
    }
    // All memory allocated on a node
    [[nodiscard]] size_t node_total_bytes() const
    {
      return node_local_bytes() + node_shmem_bytes();
    }

    // Helper functions. TODO move elsewhere and reuse?
  private:
    static std::vector<size_t>
    get_node_block_dims(const std::vector<Block_Location> &node_block_locations,
                        const std::vector<size_t> &block_dimensions)
    {
      std::vector<size_t> node_block_dims;
      node_block_dims.reserve(node_block_locations.size());
      for(auto &loc : node_block_locations)
        {
          node_block_dims.push_back(
            block_dimensions.at(loc.block_index_global));
        }
      return node_block_dims;
    }
    template <class T> static T get_max(const std::vector<T> &vec)
    {
      return *std::max_element(vec.begin(), vec.end());
    }
    template <class T> static T sum_squares(const std::vector<T> &dims)
    {
      T result(0);
      for(const auto &dim : dims)
        result += dim * dim;
      return result;
    }
  };
}