#include "Timers.hxx"
#include "sdpb_util/assert.hxx"

Scoped_Timer::Scoped_Timer(Timers &timers, const std::string &name)
    : timers(timers),
      my_timer(timers.add_and_start(name)),
      old_prefix(timers.prefix),
      new_prefix(timers.prefix + name + ".")

{
  timers.prefix = new_prefix;
}
Scoped_Timer::~Scoped_Timer()
{
  if(is_running())
    stop();
}
std::chrono::time_point<std::chrono::high_resolution_clock>
Scoped_Timer::start_time() const
{
  return my_timer.start_time;
}
void Scoped_Timer::stop()
{
  ASSERT(is_running(), "Timer '" + new_prefix + "' already stopped!");
  my_timer.stop();

  // This assertion will fail if some of the nested timers is still running:
  ASSERT(timers.prefix == new_prefix, "timers.prefix = '" + timers.prefix
                                        + "', expected: '" + new_prefix + "'");
  timers.prefix = old_prefix;
}
const Timer &Scoped_Timer::timer() const
{
  return my_timer;
}
bool Scoped_Timer::is_running() const
{
  return my_timer.is_running();
}
