//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_TIMERS_H_
#define SDPB_TIMERS_H_

#include <iostream>
#include <ostream>
#include <string>
#include <map>
#include "boost/timer/timer.hpp"

using std::map;
using std::string;
using std::ostream;
using boost::timer::cpu_timer;

// A map between strings and cpu timers
class Timers : public map<string, cpu_timer> {
 public:
  // For printing out timing information
  friend ostream& operator<<(ostream& os, const Timers& t) {
    for (map<string, cpu_timer>::const_iterator it = t.begin();
         it != t.end();
         ++it) {
      os << it->first << "\t:" << it->second.format();
    }
    return os;
  }
};

// A global Timers map for the whole program (defined in main.cpp).  A
// new timer is created by default whenver `timers' is accessed.  To
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
extern Timers timers;

#endif  // SDPB_TIMERS_H_
