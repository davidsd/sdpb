//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_TIMERS_H_
#define SDPB_TIMERS_H_

#include <iostream>
#include <ostream>
#include <map>
#include "boost/timer/timer.hpp"

using std::map;
using std::string;
using std::ostream;
using boost::timer::cpu_timer;

class Timers : public map<string,cpu_timer> {
public:
  friend ostream& operator<<(ostream& os, const Timers& t) {
    for (map<string,cpu_timer>::const_iterator it = t.begin();
         it != t.end();
         ++it) {
      os << it->first << "\t:" << it->second.format();
    }
    return os;
  }
};

extern Timers timers;

#endif  // SDPB_TIMERS_H_
