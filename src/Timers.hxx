//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Timer.hxx"

#include <El.hpp>

#include <fstream>
#include <string>
#include <list>
#include <algorithm>
#include <boost/core/noncopyable.hpp>

struct Timers : public std::list<std::pair<std::string, Timer>>
{
  bool debug = false;
  Timers(const bool &Debug) : debug(Debug) {}

  Timer &add_and_start(const std::string &name)
  {
    if(debug)
      {
        std::ifstream stat_file("/proc/self/statm");
        if(stat_file.good())
          {
            std::string stats;
            std::getline(stat_file, stats);
            El::Output(El::mpi::Rank(), " ", name, " ", stats);
          }
      }
    emplace_back(name, Timer());
    return back().second;
  }

  void write_profile(const std::string &filename) const
  {
    std::ofstream f(filename);

    f << "{" << '\n';
    for(auto it(begin()); it != end();)
      {
        f << "    {\"" << it->first << "\", " << it->second << "}";
        ++it;
        if(it != end())
          {
            f << ",";
          }
        f << '\n';
      }
    f << "}" << '\n';

    if(!f.good())
      {
        throw std::runtime_error("Error when writing to: " + filename);
      }
  }

  int64_t elapsed_milliseconds(const std::string &s) const
  {
    auto iter(std::find_if(rbegin(), rend(),
                           [&s](const std::pair<std::string, Timer> &timer) {
                             return timer.first == s;
                           }));
    if(iter == rend())
      {
        throw std::runtime_error("Could not find timing for " + s);
      }
    return iter->second.elapsed_milliseconds();
  }
};

// Simple RAII timer
// start() in constructor, stop() in destructor
struct ScopedTimer : boost::noncopyable
{
  ScopedTimer(Timers &timers, const std::string &name)
      : my_timer(timers.add_and_start(name))
  {}
  virtual ~ScopedTimer() { my_timer.stop(); }

private:
  Timer &my_timer;
};
