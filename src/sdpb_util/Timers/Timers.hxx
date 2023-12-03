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

#include <cassert>
#include <fstream>
#include <string>
#include <list>
#include <algorithm>
#include <filesystem>

struct Timers : public std::list<std::pair<std::string, Timer>>
{
  friend struct Scoped_Timer; // can change private field prefix

private:
  const bool debug = false;
  std::string prefix;

public:
  explicit Timers(bool debug);

private:
  Timer &add_and_start(const std::string &name);

public:
  void write_profile(const std::filesystem::path &path) const;

  int64_t elapsed_milliseconds(const std::string &s) const;
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
