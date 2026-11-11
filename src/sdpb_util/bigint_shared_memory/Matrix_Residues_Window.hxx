#pragma once

#include "Array_Residues_Window.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>
#include <vector>

// Vector of matrices stored in a contiguous Shared_Window_Array
// in prime-row-major order
// This is residues of Q (BLAS multiplication output is written to there)
template <class T>
struct Matrix_Residues_Window : public Array_Residues_Window<T>
{
  const size_t height;
  const size_t width;
  std::vector<El::Matrix<T>> residues;

  Matrix_Residues_Window(const El::mpi::Comm shared_memory_comm,
                         size_t num_primes, size_t height, size_t width)
      : Matrix_Residues_Window(
          Shared_Window_Array_View<T>(shared_memory_comm,
                                      num_primes * height * width),
          num_primes, height, width)
  {}
  Matrix_Residues_Window(const std::shared_ptr<Shared_Window_Array<T>> &window,
                         const size_t window_offset, size_t num_primes,
                         size_t height, size_t width)
      : Matrix_Residues_Window(
          Shared_Window_Array_View<T>(window, window_offset,
                                      num_primes * height * width),
          num_primes, height, width)
  {}
  Matrix_Residues_Window(Shared_Window_Array_View<T> window_view,
                         size_t num_primes, const size_t height,
                         const size_t width)
      : Array_Residues_Window<T>(window_view, num_primes, height * width),
        height(height),
        width(width)
  {
    ASSERT(num_primes > 0);
    ASSERT(height > 0);
    ASSERT(width > 0);
    residues.resize(num_primes);
    for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
      {
        residues.at(prime_index)
          .Attach(height, width, this->data(prime_index), height);
      }
  }
};
