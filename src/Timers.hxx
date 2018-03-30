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

// A map between strings and cpu timers
class Timers : public std::map<std::string, boost::timer::cpu_timer>
{
public:
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
    std::ofstream f(filename, std::ofstream::out | std::ofstream::trunc);

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

// A global Timers map for the whole program (defined in solve.cxx).  A
// new timer is created by default whenever `timers' is accessed.  To
// time something, simply replace
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

extern Timers timers;
