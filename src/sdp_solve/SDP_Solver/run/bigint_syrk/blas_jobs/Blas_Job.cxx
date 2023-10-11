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
}
