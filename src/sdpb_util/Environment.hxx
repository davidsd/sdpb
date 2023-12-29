#pragma once

#include <El.hpp>

struct Environment
{
  // Shared memory communicator.
  // Contains all ranks for a given node.
  El::mpi::Comm comm_shared_mem;

  Environment();
  Environment(int argc, char **argv);
  ~Environment();

  // Total number of nodes
  [[nodiscard]] int num_nodes() const;
  // Node index for current rank, 0..(num_nodes-1)
  [[nodiscard]] int node_index() const;

private:
  El::Environment env;
  int _num_nodes = -1;
  int _node_index = -1;

  void initialize();
  void finalize();
};
