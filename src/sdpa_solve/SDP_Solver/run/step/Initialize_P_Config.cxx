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
    std::vector<size_t>
    get_max_group_block_dims(const Block_Mapping &block_mapping,
                             const size_t node_index,
                             const std::vector<size_t> &block_dimensions)
    {
      std::vector<size_t> max_group_dims(block_mapping.num_groups(node_index),
                                         0);
      for(auto &loc : block_mapping.node_block_locations.at(node_index))
        {
          auto &max_group_dim = max_group_dims.at(loc.group_index);
          const auto dim = block_dimensions.at(loc.block_index_global);
          max_group_dim = std::max(max_group_dim, dim);
        }
      return max_group_dims;
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
        max_group_block_dims(get_max_group_block_dims(
          block_mapping, node_index, block_dimensions)),
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
  size_t Initialize_P_Config::G_size() const
  {
    return sum_dim_squared * max_primal_dimension_step();
  }
  size_t Initialize_P_Config::P_size() const
  {
    return sum_dim_squared * primal_dimension;
  }
  size_t Initialize_P_Config::L_normalizer_size() const
  {
    // Matrix_Normalizer stores norms = vector<BigFloat>(dim) for each block.
    // So, the total size is sum(group_size * sum(dims)) over all MPI groups on a node.
    size_t result = 0;
    for(size_t group = 0; group < block_mapping.num_groups(node_index);
        ++group)
      {
        const auto &block_map = block_mapping.mapping.at(node_index).at(group);
        for(auto &block_index : block_map.block_indices)
          {
            const auto &loc = block_mapping.block_locations.at(block_index);
            const auto dim = node_block_dims.at(loc.block_index_node);
            result += dim * block_map.num_procs;
          }
      }
    return result;
  }
  size_t Initialize_P_Config::G_normalizer_size() const
  {
    return L_normalizer_size() * max_primal_dimension_step();
  }
  size_t Initialize_P_Config::L_window_size() const
  {
    return L_size() * num_primes();
  }
  size_t Initialize_P_Config::G_window_size() const
  {
    return G_size() * num_primes();
  }
  size_t Initialize_P_Config::L_X_inv_window_bytes() const
  {
    return sizeof(double) * L_window_size();
  }
  size_t Initialize_P_Config::L_Y_window_bytes() const
  {
    return L_X_inv_window_bytes();
  }
  size_t Initialize_P_Config::G_window_bytes() const
  {
    return sizeof(double) * G_window_size();
  }
  size_t Initialize_P_Config::node_P_bytes() const
  {
    return bigfloat_bytes(precision) * P_size();
  }
  size_t Initialize_P_Config::node_L_bytes() const
  {
    return bigfloat_bytes(precision) * L_size();
  }
  size_t Initialize_P_Config::node_L_normalizer_bytes() const
  {
    return bigfloat_bytes(precision) * L_normalizer_size();
  }
  size_t Initialize_P_Config::node_G_bytes() const
  {
    return bigfloat_bytes(precision) * G_size();
  }
  size_t Initialize_P_Config::node_G_Normalizer_bytes() const
  {
    return bigfloat_bytes(precision) * G_normalizer_size();
  }
  size_t Initialize_P_Config::node_normalize_and_shift_L_bytes() const
  {
    size_t buffer_size = 0;
    for(size_t group = 0; group < block_mapping.num_groups(node_index);
        ++group)
      {
        const auto num_procs
          = block_mapping.mapping.at(node_index).at(group).num_procs;
        const auto max_dim = max_group_block_dims.at(group);
        // Max buffer size for AllReduce. Each rank allocates vector of size max_dim.
        if(num_procs > 1)
          buffer_size += max_dim * num_procs;
      }
    // Factor 3 is somewhat arbitrary bur realistic:
    // One buffer for send, one for receive, one for extra MPI stuff.
    return 3 * bigfloat_bytes(precision) * buffer_size;
  }
  size_t Initialize_P_Config::node_normalize_and_shift_G_bytes() const
  {
    return node_normalize_and_shift_L_bytes() * max_primal_dimension_step();
  }
  size_t Initialize_P_Config::node_reshape_bytes() const
  {
    // Each MPI group reshapes one block of G_p (for a single p) at a time
    // into a single column of P block.
    // We compute total memory footprint of calling Reshape() on all groups on a node.
    size_t num_elements_to_send = 0;
    for(size_t group = 0; group < block_mapping.num_groups(node_index);
        ++group)
      {
        const auto dim = max_group_block_dims.at(group);
        const auto n
          = block_mapping.mapping.at(node_index).at(group).num_procs;
        // We Reshape a (dim x dim) matrix on n ranks.
        // ~1/n of elements are local
        // Thus, we need to send ~ (n-1)/n * dim * dim elements.
        num_elements_to_send += div_ceil(dim * dim * (n - 1), n);
      }
    // Reshape() uses QueueUpdate() and ProcessQueues().
    // We need memory for sendBuf (sending local elements),
    // recvBuf (receiving local elements), and MPI buffers for AllToAll()
    // (serialized sendBuf and recvBuf data).
    // i.e. 4 buffers of roughly the same size.
    // + one for
    // sendBuf + recvBuf + serialized sendBuf + serialized recvBuf + extra MPI overhead
    // => factor of 5
    return 5 * bigfloat_bytes(precision) * num_elements_to_send;
  }

  size_t Initialize_P_Config::max_primal_dimension_step() const
  {
    return div_ceil(primal_dimension, split_factor);
  }
}