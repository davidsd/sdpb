#include "Environment.hxx"

#include "assert.hxx"

Environment::Environment() : Environment(0, nullptr) {}
Environment::Environment(int argc, char **argv) : env(argc, argv)
{
  initialize();
}
Environment::~Environment()
{
  finalize();
}
int Environment::num_nodes() const
{
  return _num_nodes;
}
int Environment::node_index() const
{
  return _node_index;
}

void Environment::initialize()
{
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
}
void Environment::finalize()
{
  _node_index = -1;
  _num_nodes = -1;
  El::mpi::Free(comm_shared_mem);
}