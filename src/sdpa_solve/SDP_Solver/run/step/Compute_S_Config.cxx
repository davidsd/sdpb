#include "Compute_S_Config.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

namespace Sdpb::Sdpa
{
  namespace
  {
    // For each MPI group on a node, return Sum(dim^2) for all blocks of this group.
    std::vector<size_t>
    total_P_block_height_per_group(const Environment &env,
                                   const Block_Info &block_info)
    {
      std::vector<size_t> result;
      const auto node_index = env.node_index();
      for(auto &group : block_info.block_mapping.mapping.at(node_index))
        {
          size_t group_height = 0;
          for(auto block_index : group.block_indices)
            {
              const auto dim = block_info.block_dimensions.at(block_index);
              group_height += dim * dim;
            }
          result.push_back(group_height);
        }
      return result;
    }

    // Find partition point of a range [begin, end) w.r.t predicate,
    // i.e. perform binary search to find p_0 such that:
    // f(p) == true for p in [begin; p_0)
    // f(p) == false for p in [p_0, end)
    // NB: fails if p_0 is out of range!
    size_t
    partition_point(size_t begin, size_t end,
                    const std::function<bool(const size_t &)> &predicate)
    {
      ASSERT(begin != end);

      const bool reverse = begin > end;
      if(reverse)
        std::swap(begin, end);

      std::vector<size_t> ps(end - begin);
      std::iota(ps.begin(), ps.end(), begin);
      if(reverse)
        std::reverse(ps.begin(), ps.end());

      const auto it = std::partition_point(ps.begin(), ps.end(), predicate);
      ASSERT(it != ps.end(), "Partition point not found!");
      return *it;
    }
  }

  // struct Compute_S_Config implementation

  size_t Compute_S_Config::node_local_bytes() const
  {
    const size_t P_bytes = initialize_P_config.P_size()
                           * bigfloat_bytes(initialize_P_config.precision);
    // initialize_P() has some temporary allocations
    // and P matrix, which persists through syrk_P
    // TODO: shall we include P_bytes in Biginit_Syrk_Config::node_local_bytes()?
    return std::max(initialize_P_config.node_local_bytes(),
                    syrk_P_config.node_local_bytes() + P_bytes);
  }
  size_t Compute_S_Config::node_shmem_bytes() const
  {
    // TODO: currently we do not reuse trmm/syrk buffers.
    // return std::max(initialize_P_config.node_shmem_bytes(),
    //                 syrk_P_config.node_shmem_bytes());
    return initialize_P_config.node_shmem_bytes()
           + syrk_P_config.node_shmem_bytes();
  }
  size_t Compute_S_Config::node_total_bytes() const
  {
    return node_local_bytes() + node_shmem_bytes();
  }

  // end struct Compute_S_Config implementation

  Compute_S_Config
  get_compute_S_config(const Environment &env, const Block_Info &block_info,
                       const size_t max_total_mem, const size_t max_shared_mem)
  {
    const auto comm = env.comm_shared_mem;
    const auto num_nodes = env.num_nodes();
    const auto precision = El::gmp::Precision();
    const size_t primal_dimension = block_info.primal_dimension;
    const auto P_group_heights
      = total_P_block_height_per_group(env, block_info);
    const size_t P_width = primal_dimension;

    const size_t max_trmm_split_factor = primal_dimension;
    const size_t max_syrk_input_split_factor
      = *std::max_element(P_group_heights.begin(), P_group_heights.end());
    const size_t max_syrk_output_split_factor = P_width;

    const auto create_init_P_cfg
      = [&](const size_t trmm_split_factor) -> Initialize_P_Config {
      return Initialize_P_Config(
        comm, precision, block_info.block_mapping, block_info.node_index,
        block_info.node_group_index(), block_info.block_dimensions,
        primal_dimension, trmm_split_factor);
    };
    const auto create_syrk_cfg
      = [&](const size_t syrk_input_split_factor,
            const size_t syrk_output_split_factor) -> Bigint_Syrk_Config {
      return Bigint_Syrk_Config(comm, precision, num_nodes, P_group_heights,
                                P_width, syrk_input_split_factor,
                                syrk_output_split_factor);
    };

    // Check that we can fit into memory with the smallest possible allocations
    // (i.e. the largest split factors)
    const Initialize_P_Config smallest_init_P_cfg
      = create_init_P_cfg(max_trmm_split_factor);
    const Bigint_Syrk_Config smallest_syrk_cfg = create_syrk_cfg(
      max_syrk_input_split_factor, max_syrk_output_split_factor);

    const Compute_S_Config smallest_cfg{smallest_init_P_cfg,
                                        smallest_syrk_cfg};

    const size_t min_total_bytes = smallest_cfg.node_total_bytes();
    const size_t min_shmem_bytes = smallest_cfg.node_shmem_bytes();
    ASSERT(min_total_bytes <= max_total_mem,
           "Not enough memory for compute_S: required at least ",
           pretty_print_bytes(max_total_mem, true),
           ", limit: ", pretty_print_bytes(min_total_bytes, true));
    ASSERT(min_shmem_bytes <= max_shared_mem,
           "Not enough shared memory for compute_S: required at least ",
           pretty_print_bytes(min_shmem_bytes, true),
           ", limit: ", pretty_print_bytes(max_shared_mem, true));

    // Now let's do binary search and find:
    // 1. minimal syrk_output_split_factor
    //   (this one is the most important since split_factor>1 is worse for performance)
    // 2. minimal syrk_input_split_factor (depends entirely on syrk_output_split_factor)
    // 3. maximal primal_dimension_step for initialize_P
    // This will allow to allocate as much memory as possible without OOM.
    // In principle, we could come up with explicit formulas for these parameters,
    // but binary search is easier and should be fast enough.
    // Every search step involves some integer arithmetics with matrix sizes,
    // and constructing Fmpz_Comb (which takes less than 1ms for prec=1024, height=1e6)

    // Minimal viable output_split_factor(syrk) = 1..primal_dimension
    size_t syrk_output_split_factor = partition_point(
      1, max_syrk_output_split_factor + 1,
      [&](const size_t output_split_factor) {
        const size_t input_split_factor = max_syrk_input_split_factor;
        const Bigint_Syrk_Config syrk_cfg
          = create_syrk_cfg(input_split_factor, output_split_factor);
        const Compute_S_Config cfg{smallest_init_P_cfg, syrk_cfg};
        return cfg.node_total_bytes() > max_total_mem
               || cfg.node_shmem_bytes() > max_shared_mem;
      });
    // NB: All nodes should have the same output split factor!
    // We need max factor to avoid OOM.
    syrk_output_split_factor = El::mpi::AllReduce(
      syrk_output_split_factor, El::mpi::MAX, El::mpi::COMM_WORLD);

    const size_t syrk_input_split_factor = partition_point(
      1, max_syrk_input_split_factor + 1,
      [&](const size_t input_split_factor) {
        const Bigint_Syrk_Config syrk_cfg
          = create_syrk_cfg(input_split_factor, syrk_output_split_factor);
        const Compute_S_Config cfg{smallest_init_P_cfg, syrk_cfg};
        return cfg.node_total_bytes() > max_total_mem
               || cfg.node_shmem_bytes() > max_shared_mem;
      });

    const Bigint_Syrk_Config syrk_P_cfg(
      comm, precision, num_nodes, P_group_heights, P_width,
      syrk_input_split_factor, syrk_output_split_factor);

    // TODO right now we are not reusing shared memory buffers.
    // Thus, Bigint_Syrk_Config will take as much shared memory as it needs,
    // and Initialize_P_Config will get the remaining amount.
    // TODO: we should either reuse buffers or e.g. give (max_shared_mem / 2) to each stage.
    const size_t trmm_split_factor = partition_point(
      1, primal_dimension + 1, [&](const size_t split_factor) {
        const Initialize_P_Config init_P_cfg = create_init_P_cfg(split_factor);
        const Compute_S_Config cfg{init_P_cfg, syrk_P_cfg};
        return cfg.node_total_bytes() > max_total_mem
               || cfg.node_shmem_bytes() > max_shared_mem;
      });
    const Initialize_P_Config initialize_P_cfg
      = create_init_P_cfg(trmm_split_factor);
    const auto result = Compute_S_Config{initialize_P_cfg, syrk_P_cfg};
    ASSERT(result.node_shmem_bytes() < max_shared_mem,
           DEBUG_STRING(result.node_shmem_bytes()),
           DEBUG_STRING(max_shared_mem));
    ASSERT(result.node_total_bytes() < max_total_mem,
           DEBUG_STRING(result.node_total_bytes()),
           DEBUG_STRING(max_total_mem));
    return result;
  }
}