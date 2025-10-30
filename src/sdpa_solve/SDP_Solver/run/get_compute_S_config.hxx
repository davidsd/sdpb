#pragma once

#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpa_solve/memory_estimates.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"
#include "step/Compute_S_Config.hxx"

namespace Sdpb::Sdpa
{
  // Estimate how many BigFloats will be allocated by SDPB on the current node,
  // (including what's already allocated, e.g. SDP)
  inline Compute_S_Config get_compute_S_config_and_print_memory(
    const Environment &env, const Block_Info &block_info, const SDP &sdp,
    const SDP_Solver &solver, const size_t max_total_memory_bytes,
    const size_t max_shared_memory_bytes, const Verbosity verbosity)
  {
    const auto &node_comm = env.comm_shared_mem;

    const auto node_reduce = [&node_comm](const size_t value) -> size_t {
      return El::mpi::Reduce(value, 0, node_comm);
    };

    // Block_Diagonal_Matrix X, Y, primal_residues, X_chol, Y_chol, dX, dY, R, Z
    const size_t X_size = node_reduce(get_matrix_size_local(solver.X));
    const size_t X_bytes = node_reduce(get_allocated_bytes(solver.X));

    // Primal_Dist_Vector x
    const size_t x_bytes = node_reduce(get_allocated_bytes(solver.x));

    // DistMatrix S ~ MxM, distributed over all nodes.
    const size_t S_size = node_reduce(get_S_size_local(sdp));
    const size_t S_bytes = node_reduce(get_S_allocated_bytes(sdp));

    // SDP struct
    // const size_t SDP_size = node_reduce(get_SDP_size_local(sdp));
    const size_t SDP_bytes = node_reduce(get_allocated_bytes(sdp));

    // SDP_Solver members
    const size_t SDP_solver_bytes = node_reduce(get_allocated_bytes(solver));

    // auto [initialize_P_bytes, trmm_bytes, trsm_bytes]
    //   = get_initialize_P_trmm_trsm_bytes(block_info, grid);
    // initialize_P_bytes = node_reduce(initialize_P_bytes);
    // trmm_bytes = node_reduce(trmm_bytes);
    // trsm_bytes = node_reduce(trsm_bytes);

    // size_t mem_total = env.node_mem_total();
    // if(mem_total == 0)
    //   mem_total = std::numeric_limits<size_t>::max();
    // mem_total -= env.initial_node_mem_used() + SDP_bytes + SDP_solver_bytes;

    const auto cfg
      = get_compute_S_config(env, block_info,
                             max_total_memory_bytes
                               - (env.initial_node_mem_used() + SDP_bytes
                                  + SDP_solver_bytes + 3 * X_bytes),
                             max_shared_memory_bytes);

    const auto compute_S_bytes = cfg.node_total_bytes();

    // run(): X_cholesky, Y_cholesky
    size_t SDP_Solver_run_bytes = 2 * X_bytes;
    // step(): dx, dX, dY, S
    SDP_Solver_run_bytes += x_bytes + 2 * X_bytes + S_bytes;
    // step():
    SDP_Solver_run_bytes += std::max(
      // Temporary allocations inside compute_S()
      compute_S_bytes,
      // Allocations after compute_S():
      // step(): minus_XY
      // compute_search_direction(): R, Z
      3 * X_bytes);

    // initial_node_mem_used() is RAM allocated at SDPB start.
    // This could be important: e.g. on 128 cores (Expanse HPC) it is ~26GB
    const size_t mem_total_bytes = env.initial_node_mem_used() + SDP_bytes
                                   + SDP_solver_bytes + SDP_Solver_run_bytes;

    if(node_comm.Rank() == 0 && verbosity >= Verbosity::debug)
      {
        // Print memory estimates

        std::vector<std::pair<std::string, size_t>> num_elements_per_category{
          {"X", X_size},
          {"S", S_size},
        };

        std::vector<std::pair<std::string, size_t>> bytes_per_category{
          {"BigFloat size", bigfloat_bytes()},
          {"Total memory estimate", mem_total_bytes},
          {"Shared memory estimate", cfg.node_shmem_bytes()},
          {"\tInitial MemUsed (at SDPB start)", env.initial_node_mem_used()},
          {"\tSDP", SDP_bytes},
          {"\tSDP_Solver", SDP_solver_bytes},
          {"\tSDP_Solver::run()", SDP_Solver_run_bytes},
        };
        if(compute_S_bytes > 3 * X_bytes)
          {
            bytes_per_category.emplace_back("\t\tcompute_S()",
                                            compute_S_bytes);
            bytes_per_category.emplace_back(
              "\t\t\tinitialize_P()",
              cfg.initialize_P_config.node_total_bytes());
            bytes_per_category.emplace_back(
              "\t\t\t\tshared memory",
              cfg.initialize_P_config.node_shmem_bytes());
            bytes_per_category.emplace_back(
              "\t\t\tsyrk_P()", cfg.syrk_P_config.node_total_bytes());
            bytes_per_category.emplace_back(
              "\t\t\t\tshared memory", cfg.syrk_P_config.node_shmem_bytes());
          }

        std::ostringstream ss;
        El::BuildStream(ss, "node=", env.node_index(),
                        " matrix sizes and memory estimates: ");

        for(const auto &[name, size] : num_elements_per_category)
          {
            El::BuildStream(ss, "\n\t#(", name, ") = ", size, " elements");
          }
        for(const auto &[name, bytes] : bytes_per_category)
          {
            El::BuildStream(ss, "\n\t", name, ": ",
                            pretty_print_bytes(bytes, true));
          }
        El::BuildStream(ss, "\n\tShared memory configuration: ",
                        "\n\t\tinitialize_P() primal dimension step: ",
                        cfg.initialize_P_config.primal_dimension_step,
                        "\n\t\tsyrk_P() input split factor: ",
                        cfg.syrk_P_config.input_split_factor,
                        "\n\t\tsyrk_P() output split factor: ",
                        cfg.syrk_P_config.output_split_factor);
        El::Output(ss.str());
      }

    return cfg;
  }

  inline Compute_S_Config get_compute_S_config_and_print_memory(
    const Environment &env, const Block_Info &block_info, const SDP &sdp,
    const SDP_Solver &solver, size_t max_shared_memory_bytes,
    const Verbosity verbosity)
  {
    size_t max_total_memory_bytes = env.node_mem_total();
    if(max_total_memory_bytes == 0)
      max_total_memory_bytes = std::numeric_limits<size_t>::max();
    // Leave some room for errors
    max_total_memory_bytes = std::ceil(0.95 * max_total_memory_bytes);

    // If user sets --maxSharedMemory limit manually, we use it.
    // Otherwise, we calculate the limit automatically.
    if(max_shared_memory_bytes == 0)
      max_shared_memory_bytes = std::numeric_limits<size_t>::max();
    return get_compute_S_config_and_print_memory(
      env, block_info, sdp, solver, max_total_memory_bytes,
      max_shared_memory_bytes, verbosity);
  }
}
