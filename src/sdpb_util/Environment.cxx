#include "Environment.hxx"

#include "Boost_Float.hxx"
#include "Proc_Meminfo.hxx"
#include "assert.hxx"

#include <csignal>

namespace
{
  bool sigterm_flag = false;

  void handle_sigterm(int signal)
  {
    ASSERT_EQUAL(signal, SIGTERM);
    sigterm_flag = true;
  }
}

Environment::Environment() : Environment(0, nullptr) {}
Environment::Environment(int argc, char **argv) : env(argc, argv)
{
  initialize();
}
Environment::~Environment()
{
  finalize();
}
void Environment::set_precision(const mp_bitcnt_t digits2)
{
  El::gmp::SetPrecision(digits2);
  // El::gmp wants base-2 bits, but boost::multiprecision wants
  // base-10 digits.
  const unsigned digits10 = El::gmp::Precision() * log(2) / log(10);
  Boost_Float::default_precision(digits10);
}
int Environment::num_nodes() const
{
  return _num_nodes;
}
int Environment::node_index() const
{
  return _node_index;
}
size_t Environment::initial_node_mem_used() const
{
  return _initial_node_mem_used;
}
// NB: This function could be made static, but it makes no sense to call it
// without constructing Environment instance and subscribing to SIGTERM.
bool Environment::sigterm_received() const
{
  return sigterm_flag;
}

void Environment::initialize()
{
  // Subscribe to SIGTERM
  std::signal(SIGTERM, handle_sigterm);

  // Create shared memory communicator (all ranks on a node)
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &comm_shared_mem.comm);

  // Determine node_index and num_nodes
  {
    int curr_node = 0;
    for(int curr_rank = 0; curr_rank < El::mpi::Size(); ++curr_rank)
      {
        // When we are on the first rank of the node, set node_id
        if(curr_rank == El::mpi::Rank() && comm_shared_mem.Rank() == 0)
          {
            _node_index = curr_node;
            curr_node++;
          }
        // All ranks on a node have the same node_index
        El::mpi::Broadcast(_node_index, 0, comm_shared_mem);
        El::mpi::Broadcast(curr_node, curr_rank, El::mpi::COMM_WORLD);
      }
    _num_nodes = curr_node;

    ASSERT(_node_index >= 0, DEBUG_STRING(_node_index));
    ASSERT(_node_index < _num_nodes, DEBUG_STRING(_node_index),
           DEBUG_STRING(_num_nodes));
  }

  // Initial MemUsed (at SDPB start)
  {
    if(comm_shared_mem.Rank() == 0)
      {
        bool res;
        const auto meminfo = Proc_Meminfo::try_read(res);
        _initial_node_mem_used = res ? meminfo.mem_used() : 0;
      }
    El::mpi::Broadcast(_initial_node_mem_used, 0, comm_shared_mem);
  }
}
void Environment::finalize()
{
  _node_index = -1;
  _num_nodes = -1;
  El::mpi::Free(comm_shared_mem);
}