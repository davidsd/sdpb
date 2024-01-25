#include <catch2/catch_amalgamated.hpp>

#include "sdpb_util/block_mapping/compute_block_grid_mapping.hxx"
#include "test_util/diff.hxx"
#include "unit_tests/util/util.hxx"

#include <unordered_set>

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

  SECTION("Random")
  {
    size_t num_nodes = GENERATE(1, 2, 5);
    size_t procs_per_node = GENERATE(1, 10, 101);
    size_t num_blocks = GENERATE(1, 2, 30, 1001);

    CAPTURE(num_nodes);
    CAPTURE(procs_per_node);
    CAPTURE(num_blocks);

    std::vector<Block_Cost> block_costs;
    {
      INFO("Generate block_costs");
      block_costs.reserve(num_blocks);
      std::default_random_engine rand_engine;
      std::uniform_int_distribution<size_t> dist(0, num_blocks * 2);
      for(size_t index = 0; index < num_blocks; ++index)
        {
          size_t cost = dist(rand_engine);
          block_costs.emplace_back(cost, index);
        }
    }

    {
      auto mapping
        = compute_block_grid_mapping(procs_per_node, num_nodes, block_costs);
      INFO("Sanity checks for compute_block_grid_mapping() result");
      // TODO add also checks for

      std::unordered_set<size_t> processed_indices;

      REQUIRE(mapping.size() == num_nodes);
      for(const auto &node_mapping : mapping)
        {
          size_t num_procs = 0;
          for(const auto &block_map : node_mapping)
            {
              REQUIRE(block_map.num_procs > 0);

              // Since each block is assigned to a single node,
              // there can be empty nodes if we don't have enough blocks.
              // Otherwise, each process should get some block.
              if(num_blocks >= num_nodes)
                REQUIRE(!block_map.block_indices.empty());

              size_t cost = 0;
              for(auto index : block_map.block_indices)
                {
                  CAPTURE(index);
                  REQUIRE(index < block_costs.size());
                  auto [it, inserted_for_a_first_time]
                    = processed_indices.insert(index);
                  REQUIRE(inserted_for_a_first_time);

                  cost += block_costs.at(index).cost;
                }
              REQUIRE(cost == block_map.cost);

              num_procs += block_map.num_procs;
            }
          INFO(
            "Mapping should contain procs_per_node processes for each node:");
          REQUIRE(num_procs == procs_per_node);
        }

      INFO("Check that all indices were added to mapping:");
      REQUIRE(processed_indices.size() == num_blocks);
    }
  }
}