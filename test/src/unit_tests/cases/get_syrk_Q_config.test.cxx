#include "catch2/catch_amalgamated.hpp"

#include "sdp_solve/SDP_Solver/run/get_syrk_Q_config.hxx"

#include <El.hpp>
#include <vector>

// get_syrk_Q_config() runs on every node with that node's memory limit and
// memory estimate. bigint_syrk_blas() reduces each submatrix Q(I,J) of the
// output split over all nodes, so the output split factor must be the same on
// all nodes even when the limits differ, e.g. with the default --maxMemory
// (a percentage of each node's MemAvailable) on nodes with different free
// memory. Before this was enforced, a node with less memory split Q's output
// window and the others didn't, and restore_and_reduce() failed with "send
// buffer has wrong size!" (16 nodes x 128 ranks, N = 11002).
TEST_CASE("get_syrk_Q_config: same output split factor on all nodes")
{
  El::mpi::Comm comm_world = El::mpi::COMM_WORLD;
  const int node_size = GENERATE(1, 2, 3);
  if(comm_world.Size() % node_size != 0)
    SKIP("number of ranks is not a multiple of node_size");
  const size_t num_nodes = comm_world.Size() / node_size;
  if(num_nodes < 2)
    SKIP("needs at least two quasi-nodes");

  // Quasi-nodes, as in calculate_Block_Matrix_square
  El::mpi::Comm node_comm;
  MPI_Comm_split(comm_world.comm, comm_world.Rank() / node_size, 0,
                 &node_comm.comm);
  const size_t node_index = comm_world.Rank() / node_size;
  REQUIRE(node_comm.Size() == node_size);

  const size_t block_width = 50; // N
  // One group per rank, heights differing from node to node
  std::vector<size_t> blocks_height_per_group(node_size,
                                              20 + 7 * node_index);
  const size_t max_group_height = blocks_height_per_group.at(0);
  const auto precision = El::gmp::Precision();
  const auto node_total_bytes = [&](const size_t input_split_factor,
                                    const size_t output_split_factor) {
    return Bigint_Syrk_Config(node_comm, precision, num_nodes,
                              blocks_height_per_group, block_width,
                              input_split_factor, output_split_factor)
      .node_total_bytes();
  };
  // Memory needed with output split factor 2 (and the largest input split)
  const size_t other_memory = 1000;
  const size_t bytes_split_2
    = node_total_bytes(max_group_height, 2) + other_memory;
  REQUIRE(bytes_split_2 < node_total_bytes(max_group_height, 1) + other_memory);

  // Only one node gets a limit that forces output split factor 2:
  // the first node (where rank 0 is), or the last one.
  const size_t tight_node = GENERATE(as<size_t>{}, 0, 1);
  const size_t tight_node_index = tight_node == 0 ? 0 : num_nodes - 1;
  DYNAMIC_SECTION("node_size=" << node_size << " num_nodes=" << num_nodes
                               << " tight node=" << tight_node_index)
  {
    const bool tight = node_index == tight_node_index;
    // Memory_Limit with 0 bytes means no limit
    const Memory_Limit max_memory(
      tight ? std::to_string(bytes_split_2) : "0", tight ? bytes_split_2 : 0);
    const Memory_Limit max_shared_memory("0", 0);

    const auto cfg = get_syrk_Q_config(node_comm, num_nodes,
                                       blocks_height_per_group, block_width,
                                       max_memory, max_shared_memory,
                                       other_memory);

    const size_t min_output_split = El::mpi::AllReduce(
      cfg.output_split_factor, El::mpi::MIN, comm_world);
    const size_t max_output_split = El::mpi::AllReduce(
      cfg.output_split_factor, El::mpi::MAX, comm_world);
    CAPTURE(node_index, cfg.output_split_factor, cfg.input_split_factor);
    // The same on all ranks, hence on all nodes
    REQUIRE(min_output_split == max_output_split);
    // The smallest factor that fits the tight node
    REQUIRE(max_output_split == 2);

    // Each node's configuration fits its own limit, with the smallest input
    // split factor allowed by the common output split factor.
    REQUIRE(cfg.node_total_bytes() + other_memory
            <= max_memory.limit_or_infinite());
    if(cfg.input_split_factor > 1)
      REQUIRE(node_total_bytes(cfg.input_split_factor - 1, 2) + other_memory
              > max_memory.limit_or_infinite());
    if(!tight)
      REQUIRE(cfg.input_split_factor == 1);
  }
}
