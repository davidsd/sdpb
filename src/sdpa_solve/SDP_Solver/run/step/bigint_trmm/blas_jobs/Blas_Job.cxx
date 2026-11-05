#include "Blas_Job.hxx"

#include "sdpb_util/assert.hxx"

namespace Sdpb::Sdpa::Trmm
{
  using Cost = Blas_Job::Cost;

  namespace
  {
    // Simple model: Cost of Trmm call is proportional to the number of multiplications.
    Cost cost_trmm(const El::LeftOrRight side, const El::Range<El::Int> &I,
                   const El::Range<El::Int> &J)
    {
      const size_t height = I.end - I.beg;
      const size_t width = J.end - J.beg;

      size_t cost;
      switch(side)
        {
        case El::LEFT: cost = height * (height + 1) / 2 * width; break;
        case El::RIGHT: cost = width * (width + 1) / 2 * height; break;
        default: LOGIC_ERROR("Unknown LeftOrRight=", side);
        }
      constexpr size_t blas_calls = 1;
      return {cost, blas_calls};
    }
  }

  Blas_Job
  Blas_Job::create_job(const El::LeftOrRight side, const El::UpperOrLower uplo,
                       const El::Orientation orientation,
                       const El::UnitOrNonUnit diag, const size_t prime_index,
                       const size_t block_index, const El::Range<El::Int> &I,
                       const El::Range<El::Int> &J)
  {
    return Blas_Job(side, uplo, orientation, diag, prime_index, block_index, I,
                    J, cost_trmm(side, I, J));
  }
  Blas_Job::Blas_Job(const El::LeftOrRight side, const El::UpperOrLower uplo,
                     const El::Orientation orientation,
                     const El::UnitOrNonUnit diag, const size_t prime_index,
                     const size_t block_index, const El::Range<El::Int> &I,
                     const El::Range<El::Int> &J, const Cost &cost)
      : side(side),
        uplo(uplo),
        orientation(orientation),
        diag(diag),
        prime_index(prime_index),
        block_index(block_index),
        I(I),
        J(J),
        cost(cost)
  {}
}