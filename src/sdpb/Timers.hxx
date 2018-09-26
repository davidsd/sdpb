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

struct Timers : public std::list<std::pair<std::string, Timer>>
{
  bool debug = false;
  Timers(const bool &Debug) : debug(Debug) {}
  
  Timer &add_and_start(const std::string &name)
  {
    if(debug)
      {
        El::Output(El::mpi::Rank(), " ", name);
      }
    emplace_back(name, Timer());
    return back().second;
  }

  void write_profile(std::string filename) const
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
