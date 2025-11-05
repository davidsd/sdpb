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
    const SDP_Solver &solver, const Solver_Parameters &parameters,
    const Verbosity verbosity)
  {
    const auto &node_comm = env.comm_shared_mem;

    const auto node_reduce = [&node_comm](const size_t value) -> size_t {
      return El::mpi::AllReduce(value, node_comm);
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

    // run(): X_cholesky, Y_cholesky
    size_t SDP_Solver_run_bytes = 2 * X_bytes;
    // step(): dx, dX, dY, S
    SDP_Solver_run_bytes += x_bytes + 2 * X_bytes + S_bytes;
    // step():
    const auto cfg = get_compute_S_config(env, block_info, parameters,
                                          SDP_bytes + SDP_solver_bytes
                                            + SDP_Solver_run_bytes);
    const auto compute_S_bytes = cfg.node_total_bytes();
    SDP_Solver_run_bytes += std::max(
      // Temporary allocations inside compute_S()
      compute_S_bytes,
      // Allocations after compute_S():
      // step(): minus_XY
      // compute_search_direction(): R, Z
      3 * X_bytes);

    // NB: this does not include initial MemUsed, which could be e.g. ~25GB on 128 cores.
    const size_t mem_required_bytes
      = SDP_bytes + SDP_solver_bytes + SDP_Solver_run_bytes;

    if(node_comm.Rank() == 0 && verbosity >= Verbosity::debug)
      {
        // Print memory estimates

        std::vector<std::pair<std::string, size_t>> num_elements_per_category{
          {"X", X_size},
          {"S", S_size},
        };

        std::vector<std::pair<std::string, size_t>> bytes_per_category{
          {"Initial MemAvailable (at SDPB start)",
           env.initial_node_mem_available()},
          {"BigFloat size", bigfloat_bytes()},
          {"Total SDPB memory estimate", mem_required_bytes},
          {"Shared memory estimate", cfg.node_shmem_bytes()},
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
                        "\n\t\tinitialize_P() trmm split factor: ",
                        cfg.initialize_P_config.split_factor,
                        "\n\t\tsyrk_P() input split factor: ",
                        cfg.syrk_P_config.input_split_factor,
                        "\n\t\tsyrk_P() output split factor: ",
                        cfg.syrk_P_config.output_split_factor);
        El::Output(ss.str());
      }

    ASSERT(mem_required_bytes <= parameters.max_memory.limit_or_infinite(),
           "Not enough memory: required ",
           pretty_print_bytes(mem_required_bytes),
           ", --maxMemory limit: ", parameters.max_memory);
    ASSERT(cfg.node_shmem_bytes()
             <= parameters.max_shared_memory.limit_or_infinite(),
           "Not enough shared memory: required ",
           pretty_print_bytes(cfg.node_shmem_bytes()),
           ", --maxSharedMemory limit: ", parameters.max_shared_memory);
    return cfg;
  }
}
