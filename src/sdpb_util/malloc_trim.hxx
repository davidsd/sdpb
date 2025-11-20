#pragma once

#include "Environment.hxx"
#include "assert.hxx"
#include "Timers/Timers.hxx"

#if defined(__GLIBC__)
#include <malloc.h>
#endif

// Force libc allocator to release deallocated memory back to OS.
// See https://man7.org/linux/man-pages/man3/malloc_trim.3.html
// NB: this is glibc-only! For example, MUSL (used in alpine-based SDPB Docker image
// ) does not have such function.
// TODO: detect other libc versions providing malloc_trim() .
inline void malloc_trim(const Environment &env, Timers &timers)
{
#if defined(__GLIBC__)
  El::mpi::Barrier(env.comm_shared_mem);

  Scoped_Timer timer(timers, "malloc_trim");
  const bool success = malloc_trim(0);
  timer.stop();

  if(!success)
    PRINT_WARNING("rank=", El::mpi::Rank(), ": malloc_trim(0) failed.");

  El::mpi::Barrier(env.comm_shared_mem);
#endif
}