#include "Blas_Job.hxx"

#include "sdpb_util/assert.hxx"

using Cost = Blas_Job::Cost;

namespace
{
  // Blas_Job calculates a submatrix of Q = P^T P.
  // Simple model for job time:
  // It is proportional to the number of output Q elements.
  // This estimate is true if we're doing naive multiplication.

  // TODO: one big BLAS call should be faster than many small ones (for the same matrix).
  // We can account for it by adding some overhead.
  // We don't know how big it is (it depends on many factors since BLAS is heavily optimized).
  // NB: infinitesimal overhead is already accounted for in LPT_scheduling(),
  // where priority queue for ranks is sorted by (cost, num_jobs).

  // calculate C_II = A_I^T A_I
  Cost cost_syrk(const El::Range<El::Int> &I)
  {
    size_t height = I.end - I.beg;

    // In syrk, we calculate only upper triangle of output matrix
    size_t elements = height * (height + 1) / 2;
    size_t blas_calls = 1;
    return {elements, blas_calls};
  }

  // calcualte C_IJ = A_I^T B_J
  Cost cost_gemm(const El::Range<El::Int> &I, const El::Range<El::Int> &J)
  {
    size_t height = I.end - I.beg;
    size_t width = J.end - J.beg;

    // In gemm, we calculate all output matrix elements
    size_t elements = height * width;
    size_t blas_calls = 1;
    return {elements, blas_calls};
  }
}

Blas_Job
Blas_Job::create_syrk_job(size_t prime_index, const El::Range<El::Int> &I)
{
  return Blas_Job(syrk, prime_index, I, I, cost_syrk(I));
}

Blas_Job
Blas_Job::create_gemm_job(size_t prime_index, const El::Range<El::Int> &I,
                          const El::Range<El::Int> &J)
{
  return Blas_Job(gemm, prime_index, I, J, cost_gemm(I, J));
}

Blas_Job::Blas_Job(Kind kind, size_t prime_index, const El::Range<El::Int> &I,
                   const El::Range<El::Int> &J, const Cost &cost)
    : kind(kind), prime_index(prime_index), I(I), J(J), cost(cost)
{
  for(const auto &range : {I, J})
    {
      ASSERT(range.beg >= 0);
      ASSERT(range.end > range.beg, DEBUG_STRING(range.beg),
             DEBUG_STRING(range.end));
    }

  // Diagonal blocks for syrk must be square!
  if(kind == syrk && I.beg == J.beg)
    ASSERT_EQUAL(I.end, J.end, DEBUG_STRING(I.beg));
}
