#include "Blas_Job.hxx"

#include "sdpb_util/assert.hxx"

using Cost = Blas_Job::Cost;

Blas_Job::Blas_Job(size_t prime_index, const El::Range<El::Int> &I,
                   const El::Range<El::Int> &J)
    : prime_index(prime_index), I(I), J(J)
{
  for(const auto &range : {I, J})
    {
      ASSERT(range.beg >= 0);
      ASSERT(range.end > range.beg);
    }

  // Diagonal blocks must be square!
  if(I.beg == J.beg)
    ASSERT_EQUAL(I.end, J.end);
}

Cost Blas_Job::cost() const
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

  size_t elements;
  // For diagonal block,
  // we calculate only upper triangle
  if(I.beg == J.beg)
    {
      ASSERT_EQUAL(height, width);
      elements = height * (height + 1) / 2;
    }
  else
    {
      // For off-diagonal blocks, we calculate all elements
      elements = height * width;
    }
  size_t blas_calls = 1;
  return {elements, blas_calls};
}

Cost::Cost() : elements(0), blas_calls(0) {}

Cost::Cost(size_t elements, size_t blas_calls)
    : elements(elements), blas_calls(blas_calls)
{}

bool Cost::operator<(const Blas_Job::Cost &rhs) const
{
  return std::tie(elements, blas_calls)
         < std::tie(rhs.elements, rhs.blas_calls);
}

Blas_Job::Cost Cost::operator+(const Blas_Job::Cost &other) const
{
  return {elements + other.elements, blas_calls + other.blas_calls};
}

Cost &Blas_Job::Cost::operator+=(const Cost &other)
{
  elements += other.elements;
  blas_calls += other.blas_calls;
  return *this;
}

std::ostream &operator<<(std::ostream &os, const Blas_Job::Cost &cost)
{
  os << "elements: " << cost.elements << " blas_calls: " << cost.blas_calls;
  return os;
}
