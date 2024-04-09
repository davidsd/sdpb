#pragma once

#include "sdpb_util/Shared_Window_Array.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>
#include <boost/noncopyable.hpp>
#include <vector>

// Vector of matrices stored in a contiguous Shared_Window_Array
// in prime-row-major order
// This is residues of Q (BLAS multiplication output is written to there)
template <class T> class Residue_Matrices_Window : boost::noncopyable
{
public:
  const size_t num_primes;
  const size_t height;
  const size_t width;
  const size_t prime_stride;
  std::vector<El::Matrix<T>> residues;

private:
  Shared_Window_Array<T> window;

public:
  Residue_Matrices_Window(El::mpi::Comm shared_memory_comm, size_t num_primes,
                          size_t height, size_t width)
      : num_primes(num_primes),
        height(height),
        width(width),
        prime_stride(height * width),
        window(shared_memory_comm, num_primes * prime_stride)
  {
    ASSERT(num_primes > 0);
    ASSERT(height > 0);
    ASSERT(width > 0);
    residues.resize(num_primes);
    for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
      {
        size_t prime_offset = prime_index * prime_stride;
        residues.at(prime_index)
          .Attach(height, width, window.data + prime_offset, height);
      }
  }
  [[nodiscard]] El::mpi::Comm Comm() const { return window.comm; }
  void Fence() { window.Fence(); }
};
