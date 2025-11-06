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
#define SCOPE(name) Scope name##_scope(#name, tracker)
#define ALLOCATION(name) Allocation name##_allocation(#name, name, tracker)
#define GROUP(name) Allocation name##_group(#name, 0, tracker)
#define GROUP_ITEM(group, name)                                               \
  Allocation name##_group(#name, name, tracker, group##_group)

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

    // Make allocations
    {
      ALLOCATION(global_buffer);
      {
        SCOPE(foo);
        ALLOCATION(foo_x);
        ALLOCATION(foo_y);
        {
          SCOPE(bar);
          GROUP(bar_xy);
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
    std::stringstream ss;
    constexpr bool also_print_exact_bytes = false;
    tracker.print(ss, also_print_exact_bytes);
    auto tracker_output = ss.str();
    // Print the whole tree
    CAPTURE(tracker_output);
    REQUIRE(tracker.is_finished());
    DIFF(tracker.peak_memory(), expected_peak_memory);
  }
}
