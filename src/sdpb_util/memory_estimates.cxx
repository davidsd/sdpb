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
      // If we don't have enough memory,
      // we want to print OOM warning for any verbosity
      if(nonshared_memory_required_per_node_bytes > mem_total_bytes
         || verbosity >= Verbosity::debug)
        {
          PRINT_WARNING(ss.str());
        }
      if(verbosity >= Verbosity::regular)
        {
          El::Output("node=", env.node_index(), ": Set --maxSharedMemory=",
                     pretty_print_bytes(max_shared_memory_bytes, true));
        }
    }
  // All ranks on a node should have the same limit
  El::mpi::Broadcast(max_shared_memory_bytes, 0, env.comm_shared_mem);

  return max_shared_memory_bytes;
}
size_t get_heap_allocated_bytes(const El::BigFloat &f)
{
  return bigfloat_bytes();
}
size_t get_heap_allocated_bytes(const El::AbstractDistMatrix<El::BigFloat> &m)
{
  return m.AllocatedMemory() * bigfloat_bytes();
}
size_t get_heap_allocated_bytes(const El::Matrix<El::BigFloat> &m)
{
  return m.MemorySize() * bigfloat_bytes();
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
