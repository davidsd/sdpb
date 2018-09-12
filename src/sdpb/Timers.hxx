//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include <boost/timer/timer.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

// Timers map.  To time something, simply replace
//
//   myComputation();
//
// with
//
//   timers["my computation"].start();
//   myComputation();
//   timers["my computation"].stop();
//
// use .resume() to accumulate time (can be used also on creation)

// A map between strings and cpu timers
class Timers
    : public std::list<std::pair<std::string, boost::timer::cpu_timer>>
{
public:
  boost::timer::cpu_timer &add_and_start(const std::string &name)
  {
    emplace_back(name, boost::timer::cpu_timer());
    return back().second;
  }

  friend std::ostream &operator<<(std::ostream &os, const Timers &timers)
  {
    unsigned long max_length = 0;
    for(auto &timer : timers)
      {
        max_length = std::max(max_length, timer.first.length());
      }
    for(auto &timer : timers)
      {
        os << std::setw(max_length) << std::left << timer.first << " :"
           << timer.second.format();
      }
    return os;
  }

  void write_profile(std::string filename)
  {
    std::ofstream f(filename);

    f << "{" << '\n';
    for(auto it(begin()); it != end();)
      {
        f << "    {\"" << it->first
          << it->second.format(10, "\", %w, %u, %s }");

        ++it;
        if(it != end())
          {
            f << ",";
          }
        f << '\n';
      }
    f << "}" << '\n';
  }
};
