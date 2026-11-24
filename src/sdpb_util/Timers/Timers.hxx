//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Timer.hxx"
#include "sdpb_util/Environment.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <El.hpp>

#include <boost/core/noncopyable.hpp>

#include <string>
#include <list>
#include <filesystem>

struct Timers
{
  friend struct Scoped_Timer; // can change private field prefix

private:
  using Named_Timers = std::list<std::pair<std::string, Timer>>;
  // We use std::list instead of std::vector:
  // each Scoped_Timer holds an iterator to its entry, and std::list keeps it
  // valid when other entries are added or erased.
  // TODO refactor timers in a way that prevents such obscure bugs.
  Named_Timers named_timers;
  std::string prefix;

  Verbosity verbosity = Verbosity::regular;
  std::string node_debug_prefix;

  bool can_read_meminfo = true;
  // Max MemUsed value
  size_t max_mem_used = 0;
  // name of the timer that had max MemUsed value
  std::string max_mem_used_name;
  // If false, a Scoped_Timer erases its entry from named_timers
  // when it is destroyed.
  bool keep_finished_timers = true;

public:
  Timers();
  Timers(const Environment &env, const Verbosity &verbosity);
  ~Timers() noexcept;

  // named_timers gets an entry for every timer started, and only
  // write_profile() reads it. Call this if write_profile() will not be called:
  // then finished timers are erased, and named_timers does not grow with the
  // number of solver iterations.
  void discard_finished_timers();

private:
  Named_Timers::iterator add_and_start(const std::string &name);

public:
  void write_profile(const std::filesystem::path &path) const;

private:
  void print_max_mem_used() const;
  void process_meminfo(const std::string &name);
};

// Simple RAII timer
// start() in constructor, stop() in destructor
// Temporarily appends name to timers.prefix
// NB: make sure that all Scoped_Timers are properly nested!
// Otherwise, prefixes will be wrong.
//
// Example:
//
// void f()
// {
//   Timers timers(false);
//   Scoped_Timer root_timer(timers, "root"); // "root"
//   Scoped_Timer foo_timer(timers, "foo"); // "root.foo"
//   foo_timer.stop();
//   Scoped_Timer bar_timer(timers, "bar"); // "root.bar"
// }
struct Scoped_Timer : boost::noncopyable
{
  Scoped_Timer(Timers &timers, const std::string &name);
  virtual ~Scoped_Timer();

  [[nodiscard]] std::chrono::time_point<std::chrono::high_resolution_clock>
  start_time() const;

  void stop();
  [[nodiscard]] const Timer &timer() const;
  [[nodiscard]] int64_t elapsed_milliseconds() const;

private:
  Timers &timers;
  // Entry in timers.named_timers
  Timers::Named_Timers::iterator my_entry;
  std::string old_prefix;
  std::string new_prefix;

  [[nodiscard]] bool is_running() const;
};
