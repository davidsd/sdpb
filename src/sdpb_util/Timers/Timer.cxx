#include "Timer.hxx"

#include <cassert>

Timer::Timer() : start_time(now()) {}
void Timer::stop()
{
  assert(is_running());
  stop_time = now();
  running = false;
}
int64_t Timer::elapsed_milliseconds() const
{
  return elapsed<std::chrono::milliseconds>();
}
int64_t Timer::elapsed_seconds() const
{
  return elapsed<std::chrono::seconds>();
}
bool Timer::is_running() const
{
  return running;
}
std::chrono::time_point<std::chrono::high_resolution_clock> Timer::now()
{
  return std::chrono::high_resolution_clock::now();
}
template <class Duration> int64_t Timer::elapsed() const
{
  auto end_time = is_running() ? now() : stop_time;
  return std::chrono::duration_cast<Duration>(end_time - start_time).count();
}
std::ostream &operator<<(std::ostream &os, const Timer &timer)
{
  os << (timer.elapsed_milliseconds() / 1000.0);
  return os;
}
