//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Timer.hxx"

#include <El.hpp>

#include <boost/core/noncopyable.hpp>

#include <string>
#include <list>
#include <filesystem>

struct Timers
{
  friend struct Scoped_Timer; // can change private field prefix

private:
  // We use std::list instead of std::vector to avoid reallocation.
  // Scoped_Timer holds reference to timer, which would be invalidated after reallocation.
  // TODO refactor timers in a way that prevents such obscure bugs.
  std::list<std::pair<std::string, Timer>> named_timers;
  const bool debug = false;
  std::string prefix;
  // Shared memory communicator, used for debug output.
  // TODO: create it somewhere near El::Environment and reuse in other places,
  El::mpi::Comm comm_shared_mem;
  bool can_read_meminfo = true;
  // Max MemUsed value
  size_t max_mem_used{};
  // name of the timer that had max MemUsed value
  std::string max_mem_used_name;

public:
  explicit Timers(bool debug);

private:
  Timer &add_and_start(const std::string &name);

public:
  void write_profile(const std::filesystem::path &path) const;

  int64_t elapsed_milliseconds(const std::string &s) const;

  void print_max_mem_used() const;

private:
  void print_meminfo(const std::string &name);
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

private:
  Timers &timers;
  Timer &my_timer;
  std::string old_prefix;
  std::string new_prefix;

  [[nodiscard]] bool is_running() const;
};
