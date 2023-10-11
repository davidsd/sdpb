#pragma once

#include <El.hpp>

// Job for calculating submatrix Q_IJ = P_I^T P_J modulo given prime
// I and J are contiguous ranges, e.g. I=[0,3), J=[6,9)
struct Blas_Job
{
  const size_t prime_index;
  const El::Range<El::Int> I;
  const El::Range<El::Int> J;

  Blas_Job(size_t prime_index, const El::Range<El::Int> &I,
           const El::Range<El::Int> &J);
};
