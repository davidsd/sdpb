#include "Initialize_P_Config.hxx"

#include "sdpb_util/assert.hxx"
#include "sdpb_util/memory_estimates.hxx"

namespace Sdpb::Sdpa
{
  // Helper functions. TODO move elsewhere and reuse?
  namespace
  {
    size_t div_ceil(const size_t a, const size_t b)
    {
      return a / b + (a % b == 0 ? 0 : 1);
    }
    std::vector<size_t>
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
    template <class T> T get_max(const std::vector<T> &vec)
    {
      return *std::max_element(vec.begin(), vec.end());
    }
    template <class T> T sum_squares(const std::vector<T> &dims)
    {
      T result(0);
      for(const auto &dim : dims)
        result += dim * dim;
      return result;
    }
  }

  Initialize_P_Config::Initialize_P_Config(const Block_Info &block_info,
                                           const size_t split_factor)
      : Initialize_P_Config(block_info.node_comm, El::gmp::Precision(),
                            block_info.block_mapping, block_info.node_index,
                            block_info.node_group_index(),
                            block_info.block_dimensions,
                            block_info.primal_dimension, split_factor)
  {}
  Initialize_P_Config::Initialize_P_Config(
    const El::mpi::Comm shared_memory_comm, const size_t precision,
    const Block_Mapping &block_mapping, const size_t node_index,
    const size_t group_index, const std::vector<size_t> &block_dimensions,
    const size_t primal_dimension, const size_t split_factor)
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
        split_factor(split_factor)
  {
    ASSERT(split_factor >= 1 && split_factor <= primal_dimension,
           DEBUG_STRING(split_factor), DEBUG_STRING(primal_dimension));
  }
  const std::vector<Block_Location> &
  Initialize_P_Config::local_block_locations() const
  {
    return block_mapping.node_group_block_locations.at(node_index)
      .at(group_index);
  }
  const std::vector<Block_Location> &
  Initialize_P_Config::node_block_locations() const
  {
    return block_mapping.node_block_locations.at(node_index);
  }
  Fmpz_Comb Initialize_P_Config::comb() const
  {
    const mp_limb_t bits = precision;
    // Max number of additions during trmm computation.
    const size_t max_additions = max_block_dim;
    // Allow negative values
    constexpr int sign = 1;
    return Fmpz_Comb(bits, bits, sign, max_additions);
  }
  size_t Initialize_P_Config::num_primes() const
  {
    return comb().num_primes;
  }
  size_t Initialize_P_Config::L_size() const
  {
    return sum_dim_squared;
  }
  size_t Initialize_P_Config::F_size() const
  {
    return sum_dim_squared * max_primal_dimension_step();
  }
  size_t Initialize_P_Config::P_size() const
  {
    return sum_dim_squared * primal_dimension;
  }
  size_t Initialize_P_Config::L_window_size() const
  {
    return L_size() * num_primes();
  }
  size_t Initialize_P_Config::F_window_size() const
  {
    return F_size() * num_primes();
  }
  size_t Initialize_P_Config::L_X_inv_window_bytes() const
  {
    return sizeof(double) * L_window_size();
  }
  size_t Initialize_P_Config::L_Y_window_bytes() const
  {
    return L_X_inv_window_bytes();
  }
  size_t Initialize_P_Config::F_window_bytes() const
  {
    return sizeof(double) * F_window_size();
  }

  size_t Initialize_P_Config::node_local_bytes() const
  {
    // Matrix_Normalizer stores norms = vector<BigFloat>(dim) for each block.
    // So, the total size is sum(group_size * sum(dims)) over all MPI groups on a node.
    size_t L_normalizer_size = 0;
    // For each block, we perform AllReduce for the vector of column norms.
    // Each MPI group is processing one block at a time.
    // So, the total buffer size <= sum(group_size * max(dim)) over all MPI groups on a node.
    size_t L_normalizer_mpi_buffer_size = 0;
    for(size_t group = 0; group < block_mapping.num_groups(node_index);
        ++group)
      {
        const auto &block_map = block_mapping.mapping.at(node_index).at(group);
        size_t max_dim = 0;
        for(auto &block_index : block_map.block_indices)
          {
            const auto &loc = block_mapping.block_locations.at(block_index);
            const auto dim = node_block_dims.at(loc.block_index_node);
            max_dim = std::max(max_dim, dim);
            L_normalizer_size += dim * block_map.num_procs;
          }
        // Max buffer size for AllReduce. Each rank allocates vector of size max_dim.
        L_normalizer_mpi_buffer_size += max_dim * block_map.num_procs;
      }

    const size_t F_normalizer_size
      = max_primal_dimension_step() * L_normalizer_size;

    // We estimate total MPI communication memory ~ 3 * buffer size.
    // The factor of 3 is somewhat arbitrary. See also comments in get_trsm_bytes().
    const size_t normalizers_communication_size
      = 3 * L_normalizer_mpi_buffer_size;

    const size_t normalizers_size = 2 * L_normalizer_size + F_normalizer_size
                                    + normalizers_communication_size;

    // L_X_inv, L_Y, F_{p_min..p_max} (also L_X_inv F_p etc.), P
    return bigfloat_bytes(precision)
           * (2 * L_size() + F_size() + P_size() + normalizers_size);
  }
  size_t Initialize_P_Config::node_shmem_bytes() const
  {
    // residues for L_X_inv, L_Y, F_p
    return sizeof(double) * (2 * L_window_size() + F_window_size());
  }
  size_t Initialize_P_Config::node_total_bytes() const
  {
    return node_local_bytes() + node_shmem_bytes();
  }
  size_t Initialize_P_Config::max_primal_dimension_step() const
  {
    return div_ceil(primal_dimension, split_factor);
  }
}