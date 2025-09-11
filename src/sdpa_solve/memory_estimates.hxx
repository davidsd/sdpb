#pragma once

#include "SDP.hxx"
#include "SDP_Solver.hxx"
#include "sdp_solve/memory_estimates.hxx"

namespace Sdpb::Sdpa
{
  inline size_t get_S_size_local(const SDP &sdp)
  {
    // #S = NxN, distributed over all nodes.
    return std::ceil(1.0 * sdp.primal_dimension() * sdp.primal_dimension()
                     / El::mpi::Size());
    // TODO: in fact, different ranks can own slightly different number of elements.
    // So this is not a precise estimate.
    // Sum of get_S_size_local() over all ranks will give upper bound to the real #(Q).

    // NB: reduce-scatter needs also ~(2 * S_size / split_factor^2 / num_nodes) per node for MPI buffers.
    // We account for it inside --maxSharedMemory limit, see BigInt_Shared_Memory_Syrk_Context() constructor.
  }
  inline size_t get_SDP_size_local(const SDP &sdp)
  {
    size_t SDP_size =
      // F_0, F_1..F_m, all have the same size
      get_matrix_size_local(sdp.sdp_block_F_0) * (sdp.primal_dimension() + 1)
      + sdp.primal_objective_c.AllocatedMemory();
    return SDP_size;
  }
  inline size_t get_heap_allocated_bytes(const SDP &sdp)
  {
    return get_heap_allocated_bytes(sdp.sdp_blocks_F)
           + get_heap_allocated_bytes(sdp.sdp_block_F_0)
           + ::get_heap_allocated_bytes(sdp.primal_objective_c);
  }
  inline size_t get_heap_allocated_bytes(const SDP_Solver &solver)
  {
    // TODO get rid of ugly global namespace prefix
    return ::get_heap_allocated_bytes(solver.x)
           + get_heap_allocated_bytes(solver.X)
           + get_heap_allocated_bytes(solver.Y)
           + ::get_heap_allocated_bytes(solver.primal_objective)
           + ::get_heap_allocated_bytes(solver.dual_objective)
           + ::get_heap_allocated_bytes(solver.duality_gap)
           + ::get_heap_allocated_bytes(solver.primal_residues)
           + ::get_heap_allocated_bytes(solver.primal_error)
           + ::get_heap_allocated_bytes(solver.dual_residues)
           + ::get_heap_allocated_bytes(solver.dual_error)
           + ::get_heap_allocated_bytes(solver.R_error);
  }
}
