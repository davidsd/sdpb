#include <cassert>
#include "Blas_Job.hxx"

Blas_Job::Blas_Job(size_t prime_index, const El::Range<El::Int> &I,
                   const El::Range<El::Int> &J)
    : prime_index(prime_index), I(I), J(J)
{
  for(const auto &range : {I, J})
    {
      assert(range.beg >= 0);
      assert(range.end > range.beg);
    }

  // Diagonal blocks must be square!
  if(I.beg == J.beg)
    assert(I.end == J.end);
}

size_t Blas_Job::cost() const
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

  size_t height = I.end - I.beg;
  size_t width = J.end - J.beg;

  // For diagonal block,
  // we calculate only upper triangle
  if(I.beg == J.beg)
    {
      assert(height == width);
      return height * (height + 1) / 2;
    }
  // For off-diagonal blocks, we calculate all elements
  return height * width;
}
