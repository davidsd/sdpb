#pragma once

#include "sdpb_util/bigint_shared_memory/blas_jobs/Blas_Job_Cost.hxx"

#include <El.hpp>

// Job for calculating submatrix Q_IJ = P_I^T P_J modulo given prime
// I and J are contiguous ranges, e.g. I=[0,3), J=[6,9)
struct Blas_Job
{
  using Cost = Blas_Job_Cost;
  enum Kind
  {
    syrk,
    gemm
  };

  Kind kind;
  size_t prime_index;
  El::Range<El::Int> I;
  El::Range<El::Int> J;
  Cost cost;

  static Blas_Job
  create_syrk_job(size_t prime_index, const El::Range<El::Int> &I);
  static Blas_Job
  create_gemm_job(size_t prime_index, const El::Range<El::Int> &I,
                  const El::Range<El::Int> &J);

private:
  Blas_Job(Kind kind, size_t prime_index, const El::Range<El::Int> &I,
           const El::Range<El::Int> &J, const Cost &cost);
};
