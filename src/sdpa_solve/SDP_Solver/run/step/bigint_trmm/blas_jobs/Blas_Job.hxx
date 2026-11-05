#pragma once

#include "sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/Blas_Job.hxx"

#include <El.hpp>

template <class TDerived> struct Abstract_Blas_Job
{
protected:
  ~Abstract_Blas_Job() = default;

public:
  Blas_Job_Cost cost;
  virtual std::vector<TDerived> split(size_t split_factor) = 0;
};

namespace Sdpb::Sdpa::Trmm
{
  // Job for calculating submatrix Q_IJ = P_I^T P_J modulo given prime
  // I and J are contiguous ranges, e.g. I=[0,3), J=[6,9)
  struct Blas_Job
  {
    using Cost = Blas_Job_Cost;

    El::LeftOrRight side;
    El::UpperOrLower uplo;
    El::Orientation orientation;
    El::UnitOrNonUnit diag;

    size_t prime_index;
    size_t block_index;
    El::Range<El::Int> I;
    El::Range<El::Int> J;
    Cost cost;

    static Blas_Job
    create_job(El::LeftOrRight side, El::UpperOrLower uplo,
               El::Orientation orientation, El::UnitOrNonUnit diag,
               size_t prime_index, size_t block_index,
               const El::Range<El::Int> &I, const El::Range<El::Int> &J);

  private:
    Blas_Job(El::LeftOrRight side, El::UpperOrLower uplo,
             El::Orientation orientation, El::UnitOrNonUnit diag,
             size_t prime_index, size_t block_index,
             const El::Range<El::Int> &I, const El::Range<El::Int> &J,
             const Cost &cost);
  };
}