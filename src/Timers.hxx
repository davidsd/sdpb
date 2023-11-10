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
      print_debug_info(name);
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

private:
  static void print_debug_info(const std::string &name)
  {
    std::ostringstream ss;
    ss << El::mpi::Rank() << " " << name << " ";
    auto prefix = ss.str();

    print_statm(prefix);

    // /proc/meminfo is the same for all processes in node,
    // so we print it only for rank 0.
    // TODO: print meminfo for a first process of each node
    // (makes sense if RAM is not distributed equally among the nodes)
    if(El::mpi::Rank() == 0)
      print_meminfo(prefix);
  }

  // TODO move to separate file

  // /proc/self/statm displays the following quantities:
  // size resident shared text lib data dt
  // Each quantity is measured in pages;
  // usually page size=4096B (you can check it by calling "getconf PAGESIZE")
  static void print_statm(const std::string &prefix)
  {
    std::ifstream stat_file("/proc/self/statm");
    if(stat_file.good())
      {
        std::string stats;
        std::getline(stat_file, stats);
        El::Output(prefix, stats);
      }
  }

  // /proc/meminfo can be different on different OSs.
  // Usually (e.g. on CentOS) it looks like
  // MemTotal:       131189996 kB
  // MemFree:        24211752 kB
  // MemAvailable:   69487008 kB
  // ...
  // We print MemAvailable (RAM available for allocation)
  // and MemUsed defined as MemUsed = MemTotal - MemAvailable.
  // MemUsed is RAM that is occupied by all processes and cannot be released
  // (i.e. it doesn't include cache)
  static void print_meminfo(const std::string &prefix)
  {
    const char *proc_meminfo_path = "/proc/meminfo";
    std::ifstream meminfo_file(proc_meminfo_path);

    if(!meminfo_file.good())
      return;

    const char *mem_total_prefix = "MemTotal:";
    const char *mem_available_prefix = "MemAvailable:";
    size_t memTotalKB = 0;
    size_t memAvailableKB = 0;
    std::string line;
    while(std::getline(meminfo_file, line))
      {
        std::istringstream iss(line);
        std::string name;
        size_t size;
        std::string kB;
        if(iss >> name >> size >> kB)
          {
            if(kB != "kB" && kB != "KB")
              {
                El::Output(proc_meminfo_path,
                           ": expected \"kB\" at the end of line: ", line);
                return;
              }
            if(name == mem_total_prefix)
              memTotalKB = size;
            else if(name == mem_available_prefix)
              memAvailableKB = size;
            if(memTotalKB > 0 && memAvailableKB > 0)
              break;
          }
        else
          {
            El::Output(proc_meminfo_path, ": cannot parse line: ", line);
            return;
          }
      }

    if(memTotalKB == 0)
      {
        El::Output(proc_meminfo_path, ": ", mem_total_prefix, " not found");
        return;
      }
    if(memAvailableKB == 0)
      {
        El::Output(proc_meminfo_path, ": ", mem_available_prefix,
                   " not found");
        return;
      }
    auto memAvailableGB = (double)memAvailableKB / 1024 / 1024;
    auto memUsedGB = (double)(memTotalKB - memAvailableKB) / 1024 / 1024;
    El::Output(prefix, "MemAvailable, GB: ", memAvailableGB);
    El::Output(prefix, "MemUsed, GB: ", memUsedGB);
  }
};
