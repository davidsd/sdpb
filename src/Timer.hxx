#pragma once

#include <cassert>
#include <chrono>
#include <iostream>

struct Timer
{
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time,
    stop_time;
  Timer() : start_time(now()) {}
  auto stop()
  {
    assert(is_running());
    stop_time = now();
    running = false;
  }

  [[nodiscard]] int64_t elapsed_milliseconds() const
  {
    return elapsed<std::chrono::milliseconds>();
  }
  [[nodiscard]] int64_t elapsed_seconds() const
  {
    return elapsed<std::chrono::seconds>();
  }

private:
  bool running = true;

  [[nodiscard]] bool is_running() const { return running; }

  static std::chrono::time_point<std::chrono::high_resolution_clock> now()
  {
    return std::chrono::high_resolution_clock::now();
  }

  template <class Duration> [[nodiscard]] int64_t elapsed() const
  {
    auto end_time = is_running() ? now() : stop_time;
    return std::chrono::duration_cast<Duration>(end_time - start_time).count();
  }
};

inline std::ostream &operator<<(std::ostream &os, const Timer &timer)
{
  os << (timer.elapsed_milliseconds() / 1000.0);
  return os;
}
