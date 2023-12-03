#pragma once

#include <chrono>
#include <iostream>

struct Timer
{
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time,
    stop_time;
  Timer();
  void stop();

  [[nodiscard]] int64_t elapsed_milliseconds() const;
  [[nodiscard]] int64_t elapsed_seconds() const;

  [[nodiscard]] bool is_running() const;

private:
  bool running = true;

  static std::chrono::time_point<std::chrono::high_resolution_clock> now();

  template <class Duration> [[nodiscard]] int64_t elapsed() const;
};

std::ostream &operator<<(std::ostream &os, const Timer &timer);
