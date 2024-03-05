#pragma once

#include <cstddef>

// This class displays data read from /proc/meminfo
//
// /proc/meminfo can be different on different OSs.
// Usually (e.g. on CentOS) it looks like
// MemTotal:       131189996 kB
// MemFree:        24211752 kB
// MemAvailable:   69487008 kB
// ...
// We store all values in bytes.
struct Proc_Meminfo
{
  // MemTotal, B
  const size_t mem_total;
  // MemAvailable, B (RAM available for allocation)
  const size_t mem_available;

  // MemUsed (bytes) is defined as MemUsed = MemTotal - MemAvailable.
  // MemUsed is RAM that is occupied by all processes and cannot be released
  // (i.e. it doesn't include cache)
  [[nodiscard]] size_t mem_used() const;

  // read from /proc/meminfo, NB: can throw exceptions
  static Proc_Meminfo read() noexcept(false);
  // read from /proc/meminfo, in case of exception catch it and set result = false
  static Proc_Meminfo
  try_read(bool &result, bool print_error_msg = false) noexcept;

private:
  Proc_Meminfo(size_t mem_total, size_t mem_available);
};
