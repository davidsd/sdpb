#pragma once

#include <chrono>
#include <iostream>

struct Timer
{
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time,
    stop_time;
  Timer() : start_time(std::chrono::high_resolution_clock::now()) {}
  auto stop() { stop_time = std::chrono::high_resolution_clock::now(); }

  int64_t elapsed_milliseconds() const
  {
    return std::chrono::duration_cast<std::chrono::milliseconds>(stop_time
                                                                 - start_time)
      .count();
  }
  int64_t elapsed_seconds() const
  {
    return std::chrono::duration_cast<std::chrono::seconds>(stop_time
                                                            - start_time)
      .count();
  }
};

inline std::ostream &operator<<(std::ostream &os, const Timer &timer)
{
  os << timer.elapsed_seconds();
  return os;
}
