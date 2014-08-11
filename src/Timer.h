#ifndef SDP_BOOTSTRAP_TIMER_H_
#define SDP_BOOTSTRAP_TIMER_H_

#include <iostream>
#include <ostream>
#include <map>
#include "boost/timer/timer.hpp"

using std::map;
using std::string;
using std::ostream;
using boost::timer::cpu_timer;

class Timer {
public:
  map<string,cpu_timer> timers;

  void start(const string &t) {
    timers[t].resume();
  }

  void stop(const string &t) {
    timers[t].stop();
  }

  friend ostream& operator<<(ostream& os, const Timer& t) {
    for (map<string,cpu_timer>::const_iterator it = t.timers.begin();
         it != t.timers.end();
         ++it) {
      os << it->first << "\t:" << it->second.format();
    }
    return os;
  }
};

extern Timer timer;

#endif  // SDP_BOOTSTRAP_TIMER_H_
