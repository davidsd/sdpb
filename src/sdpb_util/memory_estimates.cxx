#include "memory_estimates.hxx"

#include "Proc_Meminfo.hxx"
#include "assert.hxx"
#include "ostream/pretty_print_bytes.hxx"

#include <iomanip>

size_t bigfloat_bytes()
{
  return sizeof(El::BigFloat) + El::gmp::num_limbs * sizeof(mp_limb_t);
}

size_t get_max_shared_memory_bytes(
  const size_t nonshared_memory_required_per_node_bytes,
  const Environment &env, const Verbosity verbosity)
{
  El::byte can_update = false;
  // will be set only on rank=0
  size_t mem_total_bytes = 0;
  if(env.comm_shared_mem.Rank() == 0)
    {
      bool res;
      const auto meminfo
        = Proc_Meminfo::try_read(res, verbosity >= Verbosity::debug);
      mem_total_bytes = meminfo.mem_total;
      can_update = res;
    }
  El::mpi::Broadcast(can_update, 0, env.comm_shared_mem);
  if(!can_update)
    return 0;

  // Total memory required for all ranks on a node

  size_t max_shared_memory_bytes = 0;
  if(env.comm_shared_mem.Rank() == 0)
    {
      std::ostringstream ss;
      El::BuildStream(ss, "node=", env.node_index(), ": ");

      El::BuildStream(
        ss, "\n\tMemTotal: ", pretty_print_bytes(mem_total_bytes, false),
        "\n\tRequired memory estimate (excluding shared memory windows): ",
        pretty_print_bytes(nonshared_memory_required_per_node_bytes, false));

      if(nonshared_memory_required_per_node_bytes > mem_total_bytes)
        {
          // This is certainly not enough, but at least
          // we'll print sizes in BigInt_Shared_Memory_Syrk_Context
          max_shared_memory_bytes = 0.5 * mem_total_bytes;
          El::BuildStream(
            ss,
            "\n\tSDPB will set --maxSharedMemory to 50% of MemTotal, i.e. ",
            pretty_print_bytes(max_shared_memory_bytes, false));
          El::BuildStream(ss, "\n\tSDPB will probably fail with OOM. Consider "
                              "increasing number of nodes or RAM per node.");
        }
      else
        {
          // ad-hoc coefficient 0.5 to leave some free RAM
          max_shared_memory_bytes
            = 0.5
              * (mem_total_bytes - nonshared_memory_required_per_node_bytes);
          El::BuildStream(ss,
                          "\n\tTo prevent OOM, "
                          "SDPB will set --maxSharedMemory to 50% of the "
                          "remaining memory, i.e. ",
                          pretty_print_bytes(max_shared_memory_bytes, false));
          El::BuildStream(ss,
                          "\n\tIn case of OOM, consider increasing number of "
                          "nodes and/or decreasing --maxSharedMemory limit.");
        }
      if(verbosity >= Verbosity::regular)
        {
          if(env.node_index() == 0 || verbosity >= Verbosity::debug)
            PRINT_WARNING(ss.str());
        }
    }
  // All ranks on a node should have the same limit
  El::mpi::Broadcast(max_shared_memory_bytes, 0, env.comm_shared_mem);

  return max_shared_memory_bytes;
}
size_t get_matrix_size_local(const Block_Diagonal_Matrix &X)
{
  size_t X_size = 0;
  for(const auto &X_block : X.blocks)
    {
      X_size += X_block.AllocatedMemory();
    }
  return X_size;
}
size_t get_A_X_size_local(const Block_Info &block_info, const SDP &sdp)
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
size_t get_schur_complement_size_local(const Block_Info &block_info)
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
size_t get_B_size_local(const SDP &sdp)
{
  size_t B_size = 0;
  for(const auto &B_block : sdp.free_var_matrix.blocks)
    {
      B_size += B_block.AllocatedMemory();
    }
  return B_size;
}
size_t get_Q_size_local(const SDP &sdp)
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
size_t get_SDP_size_local(const SDP &sdp)
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
size_t get_heap_allocated_bytes(const El::BigFloat &f)
{
  return bigfloat_bytes();
}
size_t get_heap_allocated_bytes(const Block_Diagonal_Matrix &m)
{
  return get_heap_allocated_bytes(m.blocks);
}
size_t get_heap_allocated_bytes(const Block_Matrix &m)
{
  return get_heap_allocated_bytes(m.blocks);
}
size_t get_heap_allocated_bytes(const Block_Vector &v)
{
  return get_heap_allocated_bytes(v.blocks);
}
size_t get_heap_allocated_bytes(const El::DistMatrix<El::BigFloat> &m)
{
  return m.AllocatedMemory() * bigfloat_bytes();
}
size_t get_heap_allocated_bytes(const El::Matrix<El::BigFloat> &m)
{
  return m.MemorySize() * bigfloat_bytes();
}
size_t get_heap_allocated_bytes(const SDP &sdp)
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
size_t get_heap_allocated_bytes(const SDP_Solver &solver)
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
void print_allocation_message_per_node(const Environment &env,
                                       const std::string &name, size_t bytes)
{
  bytes = El::mpi::Reduce(bytes, 0, env.comm_shared_mem);
  if(env.comm_shared_mem.Rank() == 0)
    {
      El::Output("node=", env.node_index(), ": allocate ", name, ": ",
                 pretty_print_bytes(bytes));
    }
}
