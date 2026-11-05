#include "sdpb_util/Memory_Tracker.hxx"

#include <catch2/catch_amalgamated.hpp>

#include "test_util/diff.hxx"

#include <algorithm>
#include <unordered_set>

using Test_Util::REQUIRE_Equal::diff;

using Scope = Memory_Tracker::Scope;
using Allocation = Memory_Tracker::Allocation;

TEST_CASE("Memory_Tracker")
{
  if(El::mpi::Rank() != 0)
    return;

  INFO("Test Memory_Tracker: simulate memory allocations with different "
       "lifetimes in nested scopes and check peak memory usage.");

  {
    Memory_Tracker tracker("root");
#define SCOPE(name) Memory_Tracker::Scope name##_scope(tracker, #name)
#define ALLOCATION(name)                                                      \
  Memory_Tracker::Allocation name##_allocation(tracker, #name, name)
#define GROUP(name) Memory_Tracker::Group name##_group(tracker, #name)
#define GROUP_ITEM(group, name)                                               \
  Memory_Tracker::Allocation name##_group_item(tracker, group##_group, #name, \
                                               name)

    // Some random numbers.
    // Peak memory will occur at different locations.
    const size_t global_buffer = GENERATE(1, 100);
    const size_t foo_x = GENERATE(1, 10);
    const size_t foo_y = GENERATE(2, 20);
    const size_t bar_x = GENERATE(1, 10);
    const size_t bar_y = GENERATE(2, 20);
    const size_t bar_z = GENERATE(3, 30);
    const size_t foo_z = GENERATE(3, 30);
    const size_t baz_x = GENERATE(15, 150);

    size_t bar_xy_group_output = 0;
    // Make allocations
    {
      ALLOCATION(global_buffer);
      {
        SCOPE(foo);
        ALLOCATION(foo_x);
        ALLOCATION(foo_y);
        {
          SCOPE(bar);
          // GROUP(bar_xy);
          Memory_Tracker::Group bar_xy_group(
            tracker, "bar_xy",
            [&](const size_t peak) { bar_xy_group_output = peak; });
          GROUP_ITEM(bar_xy, bar_x);
          GROUP_ITEM(bar_xy, bar_y);
          ALLOCATION(bar_z);
        }
        ALLOCATION(foo_z);
      }
      {
        SCOPE(baz);
        ALLOCATION(baz_x);
      }
    }

    const size_t baz = baz_x;
    const size_t bar = bar_x + bar_y + bar_z;
    const size_t foo = foo_x + foo_y + std::max(bar, foo_z);

    const size_t expected_peak_memory = global_buffer + std::max(foo, baz);
    // Print the whole tree
    CAPTURE(tracker.to_string());
    REQUIRE(tracker.is_finished());
    DIFF(bar_xy_group_output, bar_x + bar_y);
    DIFF(tracker.peak_memory(), expected_peak_memory);
  }
}
