//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#pragma once

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
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
    unsigned long maxLength=0;
    for (map<string, cpu_timer>::const_iterator it = t.begin();
         it != t.end();
         ++it) {
      if(it->first.length()>maxLength) {
        maxLength=it->first.length();
      }
    }
    for (map<string, cpu_timer>::const_iterator it = t.begin();
         it != t.end();
         ++it) {
      os << std::setw(maxLength) << std::left << it->first << " :" << it->second.format(); //should be replaced with more intelligent alignment
    }
    return os;
  }

  void writeMFile (string filename)
  {
    std::ofstream f;
    f.open(filename, std::ofstream::out | std::ofstream::trunc);

    f << "{" << std::endl;

    map<string, cpu_timer>::const_iterator final = this->end();
    --final;

    for (map<string, cpu_timer>::const_iterator it = this->begin();
         it != this->end();
         ++it) {
      f << "    {\"" << it ->first << it->second.format(10, "\", %w, %u, %s }");
      if (it != final){
        f<<",";
      }
      f<<std::endl;
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

