#pragma once

#include "Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "Block_Info.hxx"
#include "SDP.hxx"
#include "SDP_Solver.hxx"
#include "sdpb_util/memory_estimates.hxx"

template<class Derived>
size_t get_matrix_size_local(const Abstract_Block_Matrix<Derived> &X)
{
  size_t X_size = 0;
  for(const auto &X_block : X.blocks)
    {
      X_size += X_block.AllocatedMemory();
    }
  return X_size;
}
inline size_t get_A_X_size_local(const Block_Info &block_info, const SDP &sdp)
{
  // A_X_inv and A_Y
  // See compute_A_X_inv()
  // Calculate on rank=0 to avoid double-counting for DistMatrices
  size_t A_X_size = 0;
  if(block_info.mpi_comm.value.Rank() == 0)
    {
      for(size_t index = 0; index < sdp.bases_blocks.size(); ++index)
        {
          const size_t block_size
            = block_info.num_points.at(block_info.block_indices.at(index / 2));
          const size_t dim
            = block_info.dimensions.at(block_info.block_indices.at(index / 2));

          for(size_t column_block = 0; column_block < dim; ++column_block)
            for(size_t row_block = 0; row_block < dim; ++row_block)
              A_X_size += block_size * block_size;
        }
    }
  return A_X_size;
}
inline size_t get_schur_complement_size_local(const Block_Info &block_info)
{
  // schur_complement + schur_complement_cholesky
  size_t schur_complement_size = 0;
  // Calculate on rank=0 to avoid double-counting for DistMatrices
  if(block_info.mpi_comm.value.Rank() == 0)
    {
      // see initialize_schur_complement_solver() and Block_Diagonal_Matrix() ctor
      const auto &block_sizes = block_info.schur_block_sizes();
      const auto &block_indices = block_info.block_indices;
      const auto num_schur_blocks = block_info.num_points.size();
      const bool scale_index = num_schur_blocks != block_sizes.size();
      for(const auto block_index : block_indices)
        {
          if(scale_index)
            {
              schur_complement_size += block_sizes.at(block_index * 2)
                                       * block_sizes.at(block_index * 2);
              schur_complement_size += block_sizes.at(block_index * 2 + 1)
                                       * block_sizes.at(block_index * 2 + 1);
            }
          else
            {
              schur_complement_size
                += block_sizes.at(block_index) * block_sizes.at(block_index);
            }
        }
    }
  return schur_complement_size;
}
inline size_t get_B_size_local(const SDP &sdp)
{
  size_t B_size = 0;
  for(const auto &B_block : sdp.free_var_matrix.blocks)
    {
      B_size += B_block.AllocatedMemory();
    }
  return B_size;
}
inline size_t get_Q_size_local(const SDP &sdp)
{
  // #Q = NxN, distributed over all nodes.
  return std::ceil(1.0 * sdp.dual_objective_b.Height()
                   * sdp.dual_objective_b.Height() / El::mpi::Size());
  // TODO: in fact, different ranks can own slightly different number of elements.
  // So this is not a precise estimate.
  // Sum of get_Q_size_local() over all ranks will give upper bound to the real #(Q).

  // NB: reduce-scatter needs also ~(2 * Q_size / split_factor^2 / num_nodes) per node for MPI buffers.
  // We account for it inside --maxSharedMemory limit, see BigInt_Shared_Memory_Syrk_Context() constructor.
}
inline size_t get_SDP_size_local(const SDP &sdp)
{
  size_t SDP_size = 0;
  for(const auto &matrix : sdp.bilinear_bases)
    SDP_size += matrix.AllocatedMemory();
  for(const auto &matrix : sdp.bases_blocks)
    SDP_size += matrix.AllocatedMemory();
  for(const auto &B_block : sdp.free_var_matrix.blocks)
    SDP_size += B_block.AllocatedMemory();
  return SDP_size;
}
template<class Derived>
size_t get_heap_allocated_bytes(const Abstract_Block_Matrix<Derived> &m)
{
  return get_heap_allocated_bytes(m.blocks);
}
inline size_t get_heap_allocated_bytes(const SDP &sdp)
{
  size_t res = get_heap_allocated_bytes(sdp.bases_blocks)
               + get_heap_allocated_bytes(sdp.bilinear_bases)
               + get_heap_allocated_bytes(sdp.dual_objective_b)
               + get_heap_allocated_bytes(sdp.free_var_matrix)
               + get_heap_allocated_bytes(sdp.objective_const)
               + get_heap_allocated_bytes(sdp.primal_objective_c);
  if(sdp.normalization.has_value())
    res += get_heap_allocated_bytes(sdp.normalization.value());
  return res;
}
inline size_t get_heap_allocated_bytes(const SDP_Solver &solver)
{
  return get_heap_allocated_bytes(solver.dual_error)
         + get_heap_allocated_bytes(solver.dual_objective)
         + get_heap_allocated_bytes(solver.dual_residues)
         + get_heap_allocated_bytes(solver.duality_gap)
         + get_heap_allocated_bytes(solver.primal_error_p)
         + get_heap_allocated_bytes(solver.primal_objective)
         + get_heap_allocated_bytes(solver.primal_residues)
         + get_heap_allocated_bytes(solver.x)
         + get_heap_allocated_bytes(solver.y)
         + get_heap_allocated_bytes(solver.primal_error_P)
         + get_heap_allocated_bytes(solver.X)
         + get_heap_allocated_bytes(solver.Y);
}
