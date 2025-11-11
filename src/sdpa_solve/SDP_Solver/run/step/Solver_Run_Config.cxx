#include "Solver_Run_Config.hxx"

#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpa_solve/memory_estimates.hxx"
#include "sdpb_util/Memory_Tracker.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

namespace Sdpb::Sdpa
{
  namespace
  {
    constexpr auto indent_string = "  ";
    constexpr bool also_print_exact_bytes = true;

    [[nodiscard]] std::string to_string(const Memory_Tracker &mem_estimates,
                                        const std::string &initial_indent = "")
    {
      return mem_estimates.to_string(initial_indent, indent_string,
                                     also_print_exact_bytes);
    }

    // For each MPI group on a node, return Sum(dim^2) for all blocks of this group.
    std::vector<size_t>
    total_P_block_height_per_group(const Block_Info &block_info)
    {
      std::vector<size_t> result;
      const auto node_index = block_info.node_index;
      for(auto &group : block_info.block_mapping.mapping.at(node_index))
        {
          size_t group_height = 0;
          for(const auto block_index : group.block_indices)
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

    Memory_Tracker
    get_memory_estimates(const Environment &env, const Solver_Run_Config &cfg,
                         size_t &sdpb_total_bytes, size_t &sdpb_shmem_bytes)
    {
      using Scope = Memory_Tracker::Scope;
      using Allocation = Memory_Tracker::Allocation;
      using Group = Memory_Tracker::Group;
      Memory_Tracker tracker("Total memory estimate");

      const auto &sizes = cfg.sizes;
      const auto &p_cfg = cfg.initialize_P_config;
      const auto &syrk_cfg = cfg.syrk_S_config;

      // Helper macros.
      // Expansion examples:
      // SCOPE(x) -> Scope x("x", tracker)
      // NAMED_SCOPE(x, "name") -> Scope x(tracker, "name")
      // FUNC_SCOPE(f) -> Scope f(tracker, "f()")
      // ALLOC(x, bytes) -> Allocation x(tracker, "x", bytes)
      // NAMED_ALLOC(x, "name", bytes) -> Allocation x(tracker, "name", bytes)
      // GROUP(x) -> Allocation x(tracker, "x", 0)
      // NAMED_GROUP(x, "name") -> Allocation x(tracker, "name", 0)
      // GROUP_ITEM(group, x, bytes) -> Allocation x(tracker, group, "x", bytes)
      // NAMED_GROUP_ITEM(group, x, "name", bytes) -> Allocation x(tracker, group, "name", bytes)
#define NAMED_SCOPE(var, ...) Scope var(tracker, __VA_ARGS__)
#define SCOPE(var, ...) NAMED_SCOPE(var, #var, __VA_ARGS__)
#define FUNC_SCOPE(func) NAMED_SCOPE(func, #func "()")
#define NAMED_ALLOC(var, name, ...) Allocation var(tracker, name, __VA_ARGS__)
#define ALLOC(var, ...) NAMED_ALLOC(var, #var, __VA_ARGS__)
      // Temporary allocation (with essentially zero lifetime).
      // Equivalent to: { Allocation tmp(tracker, name, ...); }
#define TEMP_ALLOC(name, ...)                                                 \
  std::ignore = Allocation(tracker, name, __VA_ARGS__)
#define NAMED_GROUP(var, ...) Group var(tracker, __VA_ARGS__)
#define GROUP(var, ...) NAMED_GROUP(var, #var, __VA_ARGS__)
#define NAMED_GROUP_ITEM(group, var, ...)                                     \
  Allocation var(tracker, group, __VA_ARGS__)
#define GROUP_ITEM(group, var, ...)                                           \
  NAMED_GROUP_ITEM(group, var, #var, __VA_ARGS__)
      {
        NAMED_ALLOC(init,
                    "MemUsed at SDPB start (not included in --maxMemory)",
                    env.initial_node_mem_used());
        {
          NAMED_SCOPE(sdpb, "SDPB memory allocations",
                      [&](const size_t peak) { sdpb_total_bytes = peak; });
          ALLOC(SDP, sizes.SDP_bytes);
          ALLOC(SDP_Solver, sizes.SDP_solver_bytes);
          {
            NAMED_SCOPE(run, "SDP_Solver::run()");
            ALLOC(X_cholesky, sizes.X_bytes);
            ALLOC(Y_cholesky, sizes.X_bytes);

            NAMED_GROUP(shmem, "Shared memory windows",
                        [&](const size_t peak) { sdpb_shmem_bytes = peak; });
            // initialize_P() windows
            NAMED_GROUP_ITEM(shmem, shmem_initialize_P, "initialize_P()", 0);
            GROUP_ITEM(shmem_initialize_P, L_X_inv_window,
                       p_cfg.L_X_inv_window_bytes());
            GROUP_ITEM(shmem_initialize_P, L_Y_window,
                       p_cfg.L_Y_window_bytes());
            GROUP_ITEM(shmem_initialize_P, G_window, p_cfg.G_window_bytes());
            // syrk_S() windows
            NAMED_GROUP_ITEM(shmem, shmem_syrk_S, "syrk_S()", 0);
            NAMED_GROUP_ITEM(shmem_syrk_S, syrk_input_windows,
                             "input window(s)",
                             syrk_cfg.input_windows_bytes());
            NAMED_GROUP_ITEM(shmem_syrk_S, syrk_output_window, "output window",
                             syrk_cfg.output_window_bytes());
            // step()
            {
              FUNC_SCOPE(step);
              ALLOC(dx, sizes.x_bytes);
              ALLOC(dX, sizes.X_bytes);
              ALLOC(dY, sizes.X_bytes);
              ALLOC(S, sizes.S_bytes);
              {
                FUNC_SCOPE(compute_S);
                ALLOC(P, p_cfg.node_P_bytes());
                // Temporary allocations inside initialize_P()
                {
                  FUNC_SCOPE(initialize_P);
                  ALLOC(L_X_inv, sizes.X_bytes);
                  ALLOC(L_Y, sizes.X_bytes);
                  ALLOC(L_X_inv_normalizer, p_cfg.node_L_normalizer_bytes());
                  ALLOC(L_Y_normalizer, p_cfg.node_L_normalizer_bytes());
                  TEMP_ALLOC("normalize_and_shift(L_X_inv)",
                             p_cfg.node_normalize_and_shift_L_bytes());
                  TEMP_ALLOC("normalize_and_shift(L_Y)",
                             p_cfg.node_normalize_and_shift_L_bytes());
                  {
                    NAMED_SCOPE(p_range, "loop over p index");
                    ALLOC(G, p_cfg.node_G_bytes());
                    {
                      FUNC_SCOPE(trmm);
                      ALLOC(G_normalizer, p_cfg.node_G_Normalizer_bytes());
                      TEMP_ALLOC("normalize_and_shift(G)",
                                 p_cfg.node_normalize_and_shift_G_bytes());
                    }
                    TEMP_ALLOC("reshape", p_cfg.node_reshape_bytes());
                  }
                }
                {
                  FUNC_SCOPE(syrk_S);
                  ALLOC(P_normalizer, syrk_cfg.node_P_normalizer_bytes());
                  TEMP_ALLOC("normalize_and_shift(P)",
                             syrk_cfg.node_normalize_and_shift_P_bytes());
                  TEMP_ALLOC("reduce_scatter()",
                             syrk_cfg.node_reduce_scatter_bytes());
                }
              }
              ALLOC(minus_XY, sizes.X_bytes);
              {
                FUNC_SCOPE(compute_search_direction);
                ALLOC(R, sizes.X_bytes);
                ALLOC(Z, sizes.X_bytes);
              }
            }
          }
        }
      }
      return tracker;
    }

    void print_config_and_estimates(const Environment &env,
                                    const Solver_Run_Config &cfg,
                                    const Memory_Tracker &mem_estimates,
                                    const Verbosity &verbosity)
    {
      if(verbosity < Verbosity::debug)
        return;
      if(env.comm_shared_mem.Rank() != 0)
        return;

      std::vector<std::pair<std::string, size_t>> num_elements_per_category{
        {"X", cfg.sizes.X_size},
        {"S", cfg.sizes.S_size},
      };
      std::ostringstream ss;
      size_t indent_level = 0;

#define ADD_LINE(...)                                                         \
  do                                                                          \
    {                                                                         \
      for(size_t level = 0; level < indent_level; ++level)                    \
        ss << indent_string;                                                  \
      El::BuildStream(ss, __VA_ARGS__, "\n");                                 \
  } while(false)

      ADD_LINE("node=", env.node_index(),
               " solver configuration and memory estimates:");
      ++indent_level;
      {
        ADD_LINE("Shared memory configuration:");
        ++indent_level;
        ADD_LINE("initialize_P() trmm split factor: ",
                 cfg.initialize_P_config.split_factor);
        ADD_LINE("syrk_S() input split factor: ",
                 cfg.syrk_S_config.input_split_factor);
        ADD_LINE("syrk_S() output split factor: ",
                 cfg.syrk_S_config.output_split_factor);
        --indent_level;
      }
      {
        ADD_LINE("Matrix sizes:");
        ++indent_level;
        for(const auto &[name, size] : num_elements_per_category)
          {
            ADD_LINE("#(", name, ") = ", size, " elements");
          }
        --indent_level;
      }
      ADD_LINE("BigFloat size: ",
               pretty_print_bytes(bigfloat_bytes(), also_print_exact_bytes));

      const auto initial_indent = std::string(2 * indent_level, ' ');
      mem_estimates.print(ss, initial_indent, indent_string,
                          also_print_exact_bytes);
      El::Output(ss.str());
#undef ADD_LINE
    }
  }

  Solver_Run_Config::Sizes::Sizes(const Block_Info &block_info, const SDP &sdp,
                                  const SDP_Solver &solver)
  {
    const auto node_reduce = [&block_info](const size_t value) -> size_t {
      return El::mpi::AllReduce(value, block_info.node_comm);
    };

    // Block_Diagonal_Matrix X, Y, primal_residues, X_chol, Y_chol, dX, dY, R, Z
    X_size = node_reduce(get_matrix_size_local(solver.X));
    X_bytes = node_reduce(get_allocated_bytes(solver.X));

    // Primal_Dist_Vector x
    x_bytes = node_reduce(get_allocated_bytes(solver.x));

    // DistMatrix S ~ MxM, distributed over all nodes.
    S_size = node_reduce(get_S_size_local(sdp));
    S_bytes = node_reduce(get_S_allocated_bytes(sdp));

    // SDP struct
    // const size_t SDP_size = node_reduce(get_SDP_size_local(sdp));
    SDP_bytes = node_reduce(get_allocated_bytes(sdp));

    // SDP_Solver members
    SDP_solver_bytes = node_reduce(get_allocated_bytes(solver));

    // Total height of P matrix blocks for each MPI group
    P_group_heights = total_P_block_height_per_group(block_info);
  }

  Solver_Run_Config::Solver_Run_Config(
    const Sizes &sizes, const Initialize_P_Config &initialize_P_config,
    const Bigint_Syrk_Config &syrk_S_config)
      : sizes(sizes),
        initialize_P_config(initialize_P_config),
        syrk_S_config(syrk_S_config)
  {}

  Solver_Run_Config
  Solver_Run_Config::create(const Environment &env,
                            const Block_Info &block_info, const SDP &sdp,
                            const SDP_Solver &solver,
                            const Solver_Parameters &parameters,
                            const Verbosity &verbosity)
  {
    Sizes sizes(block_info, sdp, solver);
    const auto comm = block_info.node_comm;
    const auto num_nodes = block_info.num_nodes();
    const auto precision = El::gmp::Precision();
    const size_t primal_dimension = block_info.primal_dimension;
    const auto P_group_heights = total_P_block_height_per_group(block_info);
    const size_t P_width = primal_dimension;
    const size_t max_trmm_split_factor = primal_dimension;
    const size_t max_syrk_input_split_factor = *std::max_element(
      sizes.P_group_heights.begin(), sizes.P_group_heights.end());
    const size_t max_syrk_output_split_factor = P_width;

    const size_t max_total_mem = parameters.max_memory.limit_or_infinite();
    const size_t max_shared_mem
      = parameters.max_shared_memory.limit_or_infinite();

    const auto create_init_P_cfg
      = [&](const size_t trmm_split_factor) -> Initialize_P_Config {
      return Initialize_P_Config(block_info, trmm_split_factor);
    };
    const auto create_syrk_cfg
      = [&](const size_t syrk_input_split_factor,
            const size_t syrk_output_split_factor) -> Bigint_Syrk_Config {
      return Bigint_Syrk_Config(
        comm, precision, num_nodes, sizes.P_group_heights, P_width,
        syrk_input_split_factor, syrk_output_split_factor);
    };

    // Check that we can fit into memory with the smallest possible allocations
    // (i.e. the largest split factors)
    const Initialize_P_Config smallest_init_P_cfg
      = create_init_P_cfg(max_trmm_split_factor);
    const Bigint_Syrk_Config smallest_syrk_cfg = create_syrk_cfg(
      max_syrk_input_split_factor, max_syrk_output_split_factor);

    const Solver_Run_Config smallest_cfg(sizes, smallest_init_P_cfg,
                                         smallest_syrk_cfg);
    size_t min_total_bytes = 0;
    size_t min_shmem_bytes = 0;
    const auto smallest_mem_estimates = get_memory_estimates(
      env, smallest_cfg, min_total_bytes, min_shmem_bytes);

    ASSERT(min_total_bytes <= max_total_mem,
           "Not enough memory: required at least ",
           pretty_print_bytes(min_total_bytes, true),
           ", --maxMemory limit: ", parameters.max_memory, "\n",
           to_string(smallest_mem_estimates));
    ASSERT(min_shmem_bytes <= max_shared_mem,
           "Not enough shared memory for compute_S: required at least ",
           pretty_print_bytes(min_shmem_bytes, true),
           ", --maxSharedMemory limit: ", parameters.max_shared_memory, "\n",
           to_string(smallest_mem_estimates));

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

    const auto enough_memory = [&env, &max_total_mem, &max_shared_mem, &sizes](
                                 const Initialize_P_Config &init_P_cfg,
                                 const Bigint_Syrk_Config &syrk_cfg) -> bool {
      const Solver_Run_Config cfg(sizes, init_P_cfg, syrk_cfg);
      size_t total_bytes = 0;
      size_t shmem_bytes = 0;
      get_memory_estimates(env, cfg, total_bytes, shmem_bytes);
      // Sanity check
      ASSERT(total_bytes != 0);
      ASSERT(shmem_bytes != 0);
      return total_bytes <= max_total_mem && shmem_bytes <= max_shared_mem;
    };

    // Minimal viable output_split_factor(syrk) = 1..primal_dimension
    size_t syrk_output_split_factor = partition_point(
      1, max_syrk_output_split_factor + 1,
      [&](const size_t output_split_factor) {
        const size_t input_split_factor = max_syrk_input_split_factor;
        const Bigint_Syrk_Config syrk_cfg
          = create_syrk_cfg(input_split_factor, output_split_factor);
        return !enough_memory(smallest_init_P_cfg, syrk_cfg);
      });
    // NB: All nodes should have the same output split factor!
    // We need max factor to avoid OOM.
    syrk_output_split_factor = El::mpi::AllReduce(
      syrk_output_split_factor, El::mpi::MAX, El::mpi::COMM_WORLD);

    const size_t syrk_input_split_factor
      = partition_point(1, max_syrk_input_split_factor + 1,
                        [&](const size_t input_split_factor) {
                          const Bigint_Syrk_Config syrk_cfg = create_syrk_cfg(
                            input_split_factor, syrk_output_split_factor);
                          return !enough_memory(smallest_init_P_cfg, syrk_cfg);
                        });

    const Bigint_Syrk_Config syrk_S_cfg(
      comm, precision, num_nodes, P_group_heights, P_width,
      syrk_input_split_factor, syrk_output_split_factor);

    // TODO right now we are not reusing shared memory buffers.
    // Thus, Bigint_Syrk_Config will take as much shared memory as it needs,
    // and Initialize_P_Config will get the remaining amount.
    // TODO: we should either reuse buffers or e.g. give (max_shared_mem / 2) to each stage.
    const size_t trmm_split_factor = partition_point(
      1, primal_dimension + 1, [&](const size_t split_factor) {
        const Initialize_P_Config init_P_cfg = create_init_P_cfg(split_factor);
        return !enough_memory(init_P_cfg, syrk_S_cfg);
      });
    const Initialize_P_Config initialize_P_cfg
      = create_init_P_cfg(trmm_split_factor);
    const auto result = Solver_Run_Config(sizes, initialize_P_cfg, syrk_S_cfg);

    {
      size_t node_total_sdpb_bytes = 0;
      size_t node_shmem_sdpb_bytes = 0;
      const auto mem_estimates = get_memory_estimates(
        env, result, node_total_sdpb_bytes, node_shmem_sdpb_bytes);

      print_config_and_estimates(env, result, mem_estimates, verbosity);

      ASSERT(node_total_sdpb_bytes < max_total_mem,
             DEBUG_STRING(node_total_sdpb_bytes),
             DEBUG_STRING(parameters.max_memory),
             DEBUG_STRING(block_info.node_index), "\n",
             to_string(mem_estimates));
      ASSERT(node_shmem_sdpb_bytes < max_shared_mem,
             DEBUG_STRING(node_shmem_sdpb_bytes),
             DEBUG_STRING(parameters.max_shared_memory),
             DEBUG_STRING(block_info.node_index), "\n",
             to_string(mem_estimates));
    }
    return result;
  }
}