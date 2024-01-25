#include <catch2/catch_amalgamated.hpp>

#include "sdpb_util/block_mapping/compute_block_grid_mapping.hxx"
#include "test_util/diff.hxx"
#include "unit_tests/util/util.hxx"

using Test_Util::REQUIRE_Equal::diff;

void diff(const Block_Map &a, const Block_Map &b)
{
  DIFF(a.num_procs, b.num_procs);
  DIFF(a.cost, b.cost);
  DIFF(a.block_indices, b.block_indices);
}

TEST_CASE("compute_block_grid_mapping")
{
  if(El::mpi::Rank() != 0)
    return;

  SECTION("Simple unsorted costs")
  {
    // This simple test failed because
    // compute_block_grid_mapping() did not sort block costs.
    // TODO: add more tests
    const size_t procs_per_node = 3;
    const size_t num_nodes = 1;
    std::vector<Block_Cost> block_costs = {{1, 0}, {2, 1}};
    auto mapping
      = compute_block_grid_mapping(procs_per_node, num_nodes, block_costs);

    // TODO in principle, order of Block_Map's in expected result is not fixed.
    // Here we choose the one that is generated
    // by current compute_block_grid_mapping() implementation.
    const decltype(mapping) expected_mapping
      = {{Block_Map(2, 2, {1}), Block_Map(1, 1, {0})}};
    DIFF(mapping, expected_mapping);

    {
      std::sort(block_costs.rbegin(), block_costs.rend());
      mapping
        = compute_block_grid_mapping(procs_per_node, num_nodes, block_costs);
      DIFF(mapping, expected_mapping);
    }
  }
}