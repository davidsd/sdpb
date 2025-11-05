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

  // Set precision for El::BigFloat (GMP) and Boost_Float (MPFR)
  static void set_precision(mp_bitcnt_t digits2);

  // Total number of nodes
  [[nodiscard]] int num_nodes() const;
  // Node index for current rank, 0..(num_nodes-1)
  [[nodiscard]] int node_index() const;

  // Memory available on a node at the initialization (in bytes)
  [[nodiscard]] size_t initial_node_mem_available() const;
  // Memory used on a node at the initialization (in bytes)
  [[nodiscard]] size_t initial_node_mem_used() const;

  // Total Memory on a node(in bytes)
  [[nodiscard]] size_t node_mem_total() const;

  [[nodiscard]] bool sigterm_received() const;

private:
  El::Environment env;
  size_t _initial_node_mem_used = 0;
  size_t _node_mem_total = 0;
  int _num_nodes = -1;
  int _node_index = -1;

  void initialize();
  void finalize();
};
