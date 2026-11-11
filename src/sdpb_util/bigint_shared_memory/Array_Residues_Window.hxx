#pragma once

#include "sdpb_util/Shared_Window_Array.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

// Residues of an array of e.g. BigInts, stored as:
// array[i] % prime[prime_index] -> window[prime_index * prime_stride + i]
template <class T> struct Array_Residues_Window
{
  const size_t num_primes;
  const size_t prime_stride;
  Shared_Window_Array_View<T> window_view;

  Array_Residues_Window(Shared_Window_Array_View<T> window_view,
                        const size_t num_primes, const size_t prime_stride)
      : num_primes(num_primes),
        prime_stride(prime_stride),
        window_view(window_view)
  {
    // TODO: check for strict equality? This would require adjustments in calling code.
    // ASSERT_EQUAL(window_view.size(), num_primes * prime_stride,
    //              DEBUG_STRING(num_primes), DEBUG_STRING(prime_stride));
    ASSERT(window_view.size() >= num_primes * prime_stride,
           DEBUG_STRING(window_view.size()), DEBUG_STRING(num_primes),
           DEBUG_STRING(prime_stride));
  }

  size_t prime_offset(const size_t prime_index) const
  {
    return prime_index * prime_stride;
  }

  T *data(const size_t prime_index)
  {
    return window_view.data() + prime_offset(prime_index);
  }
  [[nodiscard]] El::mpi::Comm Comm() const { return window_view.Comm(); }
  void Fence() { window_view.Fence(); }
  [[nodiscard]] size_t size() const { return window_view.size(); }
};
