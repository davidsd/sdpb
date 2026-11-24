#include <unistd.h>
#include "sdpb_util/Timers/Timers.hxx"
#include "unit_tests/util/util.hxx"

#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>

// What a profile holds, with and without Timers::discard_finished_timers().
// A temporary file per rank, removed at the end.
namespace
{
  std::string profile_of(const Timers &timers)
  {
    const auto path = std::filesystem::temp_directory_path()
                      / ("Timers.test." + std::to_string(El::mpi::Rank()) + "."
                         + std::to_string(::getpid()));
    timers.write_profile(path);
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    std::filesystem::remove(path);
    return ss.str();
  }
  bool has(const std::string &profile, const std::string &name)
  {
    return profile.find("\"" + name + "\"") != std::string::npos;
  }
}

TEST_CASE("Timers")
{
  SECTION("by default, every finished timer stays in the profile")
  {
    Timers timers;
    {
      Scoped_Timer outer(timers, "outer");
      Scoped_Timer inner(timers, "inner");
    }
    const auto profile = profile_of(timers);
    REQUIRE(has(profile, "outer"));
    REQUIRE(has(profile, "outer.inner"));
  }
  SECTION("discard_finished_timers: finished timers go, running ones stay")
  {
    Timers timers;
    timers.discard_finished_timers();
    Scoped_Timer outer(timers, "outer");
    for(int i = 0; i < 3; ++i)
      {
        Scoped_Timer iter(timers, "iter_" + std::to_string(i));
        Scoped_Timer step(timers, "step");
        // A timer stopped before it is destroyed goes too
        step.stop();
      }
    const auto profile = profile_of(timers);
    REQUIRE(has(profile, "outer"));
    REQUIRE_FALSE(has(profile, "outer.iter_0"));
    REQUIRE_FALSE(has(profile, "outer.iter_2.step"));
  }
  SECTION("discard_finished_timers: a timer destroyed before a later sibling")
  {
    // Scoped_Timers need not be destroyed in the reverse order of creation:
    // each Scoped_Timer must erase its own entry, not the last one.
    Timers timers;
    timers.discard_finished_timers();
    std::optional<Scoped_Timer> first, second;
    first.emplace(timers, "first");
    first->stop();
    second.emplace(timers, "second");
    REQUIRE(has(profile_of(timers), "first"));

    first.reset();
    auto profile = profile_of(timers);
    REQUIRE_FALSE(has(profile, "first"));
    REQUIRE(has(profile, "second"));

    second.reset();
    REQUIRE_FALSE(has(profile_of(timers), "second"));
  }
}
