#pragma once

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

// fn(i) for 0 <= i < n on OMP_NUM_THREADS threads (one if it is unset), each
// taking the next i, so that blocks of very different sizes balance. Not OpenMP:
// Elemental's headers change with _OPENMP, and they must read alike in all of
// SDPB. fn must not touch MPI. The first exception is rethrown on this thread.
template <class Fn> void for_each_threaded(const size_t n, const Fn &fn)
{
  const char *s = std::getenv("OMP_NUM_THREADS");
  const size_t nt = std::min<size_t>(n, s ? std::max(1, std::atoi(s)) : 1);
  std::atomic<size_t> next(0);
  std::exception_ptr error;
  std::mutex error_mutex;
  auto work = [&]() {
    for(size_t i = next++; i < n; i = next++)
      {
        try
          {
            fn(i);
          }
        catch(...)
          {
            std::lock_guard<std::mutex> lock(error_mutex);
            if(!error)
              error = std::current_exception();
            next = n;
          }
      }
  };
  std::vector<std::thread> pool;
  for(size_t t = 1; t < nt; ++t)
    pool.emplace_back(work);
  work();
  for(auto &th : pool)
    th.join();
  if(error)
    std::rethrow_exception(error);
}
