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
  // For printing out timing information
  friend std::ostream &operator<<(std::ostream &os, const Timers &t)
  {
    unsigned long maxLength = 0;
    for(std::map<std::string, boost::timer::cpu_timer>::const_iterator it
        = t.begin();
        it != t.end(); ++it)
      {
        if(it->first.length() > maxLength)
          {
            maxLength = it->first.length();
          }
      }
    for(std::map<std::string, boost::timer::cpu_timer>::const_iterator it
        = t.begin();
        it != t.end(); ++it)
      {
        os << std::setw(maxLength) << std::left << it->first << " :"
           << it->second.format(); // should be replaced with more intelligent
                                   // alignment
      }
    return os;
  }

  void writeMFile(std::string filename)
  {
    std::ofstream f;
    f.open(filename, std::ofstream::out | std::ofstream::trunc);

    f << "{" << std::endl;

    std::map<std::string, boost::timer::cpu_timer>::const_iterator final
      = this->end();
    --final;

    for(std::map<std::string, boost::timer::cpu_timer>::const_iterator it
        = this->begin();
        it != this->end(); ++it)
      {
        f << "    {\"" << it->first
          << it->second.format(10, "\", %w, %u, %s }");
        if(it != final)
          {
            f << ",";
          }
        f << std::endl;
      }

    f << "}" << std::endl;

    f.close();
  }
};

// A global Timers map for the whole program (defined in main.cpp).  A
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
