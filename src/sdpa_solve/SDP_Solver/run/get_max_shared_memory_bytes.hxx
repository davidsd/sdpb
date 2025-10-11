#pragma once

#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpa_solve/memory_estimates.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

namespace Sdpb::Sdpa
{
  // Memory allocated in initialize_P(), see in step/compute_S.cxx
  // Includes:
  // - Size of matrix P (Block_Matrix part on current node)
  // - Memory usage for processing the largest block in a group:
  //   - F_p L_Y block
  //   - G_p = L_X_inv F_p L_Y.
  //   - Extra memory usage inside El::Trsm() call.
  // NB: returns all memory for MPI group from group_rank = 0
  // Other ranks return 0, so that mpi::Reduce() returns correct total memory.
  inline size_t
  get_initialize_P_bytes(const Block_Info &block_info, const El::Grid &grid)
  {
    const auto &group_comm = block_info.mpi_comm.value;
    const size_t primal_dimension = block_info.primal_dimension;

    if(group_comm.Rank() != 0)
      return 0;

    // Height of P submatrix stored on current MPI group
    size_t P_group_height = 0;
    const size_t P_width = primal_dimension;

    // Choose the largest block in the group
    size_t max_block_dim = 0;
    for(const auto &block_index : block_info.block_indices)
      {
        const size_t dim = block_info.block_dimensions.at(block_index);
        P_group_height += dim * dim;
        max_block_dim = std::max(max_block_dim, dim);
      }

    const size_t P_size = P_group_height * P_width;
    const size_t P_bytes = bigfloat_bytes() * P_size;

    const size_t FY_block_vertical_bytes
      = bigfloat_bytes() * max_block_dim * max_block_dim * primal_dimension;
    const size_t G_block_horizontal_bytes = FY_block_vertical_bytes;

    const size_t trsm_bytes
      = get_trsm_bytes(max_block_dim, max_block_dim * primal_dimension,
                       grid.Height(), grid.Width());
    return P_bytes + FY_block_vertical_bytes + G_block_horizontal_bytes
           + trsm_bytes;
  }

  // Estimate how many BigFloats will be allocated by SDPB on the current node,
  // (including what's already allocated, e.g. SDP)
  inline size_t get_required_nonshared_memory_per_node_bytes(
    const Environment &env, const Block_Info &block_info, const SDP &sdp,
    const El::Grid &grid, const SDP_Solver &solver, const Verbosity verbosity)
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

    const size_t initialize_P_bytes
      = node_reduce(get_initialize_P_bytes(block_info, grid));

    // run(): X_cholesky, Y_cholesky
    size_t SDP_Solver_run_bytes = 2 * X_bytes;
    // step(): dx, dX, dY, S
    SDP_Solver_run_bytes += x_bytes + 2 * X_bytes + S_bytes;
    // step():
    SDP_Solver_run_bytes += std::max(
      // Temporary allocations inside compute_S()
      initialize_P_bytes,
      // Allocations after compute_S():
      // step(): minus_XY
      // compute_search_direction(): R, Z
      3 * X_bytes);

    // We will use only result on rank=0
    if(node_comm.Rank() != 0)
      return 0;

    // Calculate mem_required_size

    // initial_node_mem_used() is RAM allocated at SDPB start.
    // This could be important: e.g. on 128 cores (Expanse HPC) it is ~26GB
    const size_t mem_required_bytes = env.initial_node_mem_used() + SDP_bytes
                                      + SDP_solver_bytes
                                      + SDP_Solver_run_bytes;

    if(verbosity >= Verbosity::debug)
      {
        // Print memory estimates

        std::vector<std::pair<std::string, size_t>> num_elements_per_category{
          {"X", X_size},
          {"S", S_size},
        };

        std::vector<std::pair<std::string, size_t>> bytes_per_category{
          {"BigFloat size", bigfloat_bytes()},
          {"Total non-shared memory estimate", mem_required_bytes},
          {"\tInitial MemUsed (at SDPB start)", env.initial_node_mem_used()},
          {"\tSDP", SDP_bytes},
          {"\tSDP_Solver", SDP_solver_bytes},
          {"\tSDP_Solver::run()", SDP_Solver_run_bytes},
        };
        if(initialize_P_bytes > 3 * X_bytes)
          {
            bytes_per_category.emplace_back("\t\tinitialize_P()",
                                            initialize_P_bytes);
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
        El::Output(ss.str());
      }

    return mem_required_bytes;
  }

  inline size_t
  get_max_shared_memory_bytes(const size_t default_max_shared_memory_bytes,
                              const Environment &env,
                              const Block_Info &block_info, const SDP &sdp,
                              const El::Grid &grid, const SDP_Solver &solver,
                              const Verbosity verbosity)
  {
    // If user sets --maxSharedMemory limit manually, we use it.
    // Otherwise, we calculate the limit automatically.
    if(default_max_shared_memory_bytes != 0)
      return default_max_shared_memory_bytes;
    const size_t nonshared_memory_required_per_node_bytes
      = get_required_nonshared_memory_per_node_bytes(env, block_info, sdp,
                                                     grid, solver, verbosity);
    return get_max_shared_memory_bytes(
      nonshared_memory_required_per_node_bytes, env, verbosity);
  }
}
