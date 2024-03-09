#pragma once

#include <El.hpp>

// Job for calculating submatrix Q_IJ = P_I^T P_J modulo given prime
// I and J are contiguous ranges, e.g. I=[0,3), J=[6,9)
struct Blas_Job
{
  enum Kind
  {
    syrk,
    gemm
  };

  struct Cost
  {
    size_t elements;
    size_t blas_calls;
    Cost();
    Cost(size_t elements, size_t blas_calls);
    bool operator<(const Cost &rhs) const;
    Cost operator+(Cost const &other) const;
    Cost &operator+=(const Cost &other);
    friend std::ostream &operator<<(std::ostream &os, const Cost &cost);
  };

  const Kind kind;
  const size_t prime_index;
  const El::Range<El::Int> I;
  const El::Range<El::Int> J;
  const Cost cost;

  static Blas_Job
  create_syrk_job(size_t prime_index, const El::Range<El::Int> &I);
  static Blas_Job
  create_gemm_job(size_t prime_index, const El::Range<El::Int> &I,
                  const El::Range<El::Int> &J);

private:
  Blas_Job(Kind kind, size_t prime_index, const El::Range<El::Int> &I,
           const El::Range<El::Int> &J, const Cost &cost);
};
