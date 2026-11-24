#include "Timers.hxx"
#include "sdpb_util/assert.hxx"

Scoped_Timer::Scoped_Timer(Timers &timers, const std::string &name)
    : timers(timers),
      my_entry(timers.add_and_start(name)),
      old_prefix(timers.prefix),
      new_prefix(timers.prefix + name + ".")

{
  timers.prefix = new_prefix;
}
Scoped_Timer::~Scoped_Timer()
{
  if(is_running())
    stop();
  if(!timers.keep_finished_timers)
    timers.named_timers.erase(my_entry);
}
std::chrono::time_point<std::chrono::high_resolution_clock>
Scoped_Timer::start_time() const
{
  return my_entry->second.start_time;
}
void Scoped_Timer::stop()
{
  ASSERT(is_running(), "Timer '" + new_prefix + "' already stopped!");
  my_entry->second.stop();

  // This assertion will fail if some of the nested timers is still running:
  ASSERT_EQUAL(timers.prefix, new_prefix);
  timers.prefix = old_prefix;
}
const Timer &Scoped_Timer::timer() const
{
  return my_entry->second;
}
bool Scoped_Timer::is_running() const
{
  return my_entry->second.is_running();
}
int64_t Scoped_Timer::elapsed_milliseconds() const
{
  return my_entry->second.elapsed_milliseconds();
}
