#pragma once

#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

inline BigInt_Shared_Memory_Syrk_Context
initialize_bigint_syrk_context(const Block_Info &block_info, const SDP &sdp,
                               const size_t max_shared_memory_bytes,
                               const Verbosity verbosity)
{
  const auto group_index = block_info.node_group_index();
  const auto precision = El::gmp::Precision();
  const auto num_nodes = block_info.num_nodes();
  const auto &node_comm = block_info.node_comm;

  const auto block_width = sdp.dual_objective_b.Height(); // = N

  std::vector<size_t> blocks_height_per_group;
  size_t max_group_height = 0;
  for(const auto &group :
      block_info.block_mapping.mapping.at(block_info.node_index))
    {
      size_t height = 0;
      for(const auto &block_index : group.block_indices)
        {
          // P' = height of the current block of sdp.free_var_matrix (aka B)
          height += block_info.get_schur_block_size(block_index);
        }
      blocks_height_per_group.push_back(height);
      max_group_height = std::max(max_group_height, height);
    }

  const size_t max_input_split_factor = max_group_height;
  const size_t max_output_split_factor = block_width;

  const Bigint_Syrk_Config smallest_cfg(
    node_comm, precision, num_nodes, blocks_height_per_group, block_width,
    max_input_split_factor, max_output_split_factor);
  const Bigint_Syrk_Config largest_cfg(node_comm, precision, num_nodes,
                                       blocks_height_per_group, block_width, 1,
                                       1);

  const size_t min_shmem_bytes = smallest_cfg.node_shmem_bytes();
  ASSERT(min_shmem_bytes <= max_shared_memory_bytes,
         "Not enough shared memory for compute_Q(): required at least ",
         pretty_print_bytes(min_shmem_bytes, true),
         ", limit: ", pretty_print_bytes(max_shared_memory_bytes, true));

  // Find partition point of a range [begin, end) w.r.t predicate,
  // i.e. perform binary search to find p_0 such that:
  // f(p) == true for p in [begin; p_0)
  // f(p) == false for p in [p_0, end)
  // NB: fails if p_0 is out of range!
  // TODO deduplicate with Compute_S_Config.cxx
  constexpr auto partition_point_unsafe
    = [](size_t begin, size_t end,
         const std::function<bool(const size_t &)> &predicate) -> size_t {
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
  };

  const size_t output_split_factor = partition_point_unsafe(
    1, max_output_split_factor + 1, [&](const size_t split_factor) {
      const size_t input_split_factor = max_input_split_factor;
      const Bigint_Syrk_Config cfg(node_comm, precision, num_nodes,
                                   blocks_height_per_group, block_width,
                                   input_split_factor, split_factor);
      return cfg.node_shmem_bytes() > max_shared_memory_bytes;
    });

  const size_t input_split_factor = partition_point_unsafe(
    1, max_input_split_factor + 1, [&](const size_t split_factor) {
      const Bigint_Syrk_Config cfg(node_comm, precision, num_nodes,
                                   blocks_height_per_group, block_width,
                                   split_factor, output_split_factor);
      return cfg.node_shmem_bytes() > max_shared_memory_bytes;
    });

  const Bigint_Syrk_Config cfg(node_comm, precision, num_nodes,
                               blocks_height_per_group, block_width,
                               input_split_factor, output_split_factor);

  // Print warnings for large split factors.
  if(node_comm.Rank() == 0)
    {
      const bool print_output_warning = cfg.output_split_factor > 1;
      const bool print_input_warning = cfg.input_split_factor > 10;
      if(print_output_warning || print_input_warning)
        {
          std::ostringstream ss;
          El::BuildStream(ss, "rank=", El::mpi::Rank(),
                          ": BigInt_Shared_Memory_Syrk_Context:\n");
          if(print_output_warning)
            {
              El::BuildStream(ss, "\tOutput window is split by a factor of ",
                              cfg.output_split_factor,
                              ", which may affect performance.\n");
            }
          if(print_input_warning)
            {
              El::BuildStream(
                ss, "\tInput window is split by a large factor of ",
                cfg.input_split_factor, ", which may affect performance.\n");
            }
          El::BuildStream(
            ss,
            "\tConsider increasing available shared memory per node "
            "(--maxSharedMemory option).",
            "\n\tShared memory limit: ",
            pretty_print_bytes(max_shared_memory_bytes, true),
            "\n\tOutput window:", "\n\t\tactual size after splitting: ",
            pretty_print_bytes(cfg.output_window_size() * sizeof(double),
                               true),
            "\n\t\toptimal size without splitting: ",
            pretty_print_bytes(
              largest_cfg.output_window_size() * sizeof(double), true),
            "\n\tInput window:"
            "\n\t\tactual size after splitting: ",
            pretty_print_bytes(cfg.input_window_size() * sizeof(double), true),
            "\n\t\toptimal size without splitting: ",
            pretty_print_bytes(
              largest_cfg.input_window_size() * sizeof(double), true));
          PRINT_WARNING(ss.str());
        }
    }

  return BigInt_Shared_Memory_Syrk_Context(
    cfg, group_index, block_info.block_indices, verbosity);
}
