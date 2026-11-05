#include "sdpb_util/bigint_shared_memory/blas_jobs/create_blas_job_schedule.hxx"

#include "create_blas_job_schedule.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/split_range.hxx"
#include "sdpb_util/bigint_shared_memory/blas_jobs/Blas_Job_Schedule.hxx"

namespace Sdpb::Sdpa::Trmm
{
  // Jobs for X =: L X (or X := X L) modulo each prime
  // L is Block_Diagonal_Matrix
  // X is vector<Block_Diagonal_Matrix>, vector length = vector_size
  // L and each element of X have the same block structure defined by block_dims
  //
  // If side == El::LEFT
  //   X := L X
  //   All elements of vector X are stacked horizontally,
  //   i.e. each block of X is now wide, ~ M_i x (M_i*vector_size)
  // If side == El::RIGHT:
  //   X := X L
  //   All elements of vector X are stacked vertically,
  //   i.e. each block of X is now tall, ~ (M_i * vector_size) x M_i
  // Each block of L_i is triangular, ~ M_i x M_i.
  Blas_Job_Schedule<Blas_Job>
  create_blas_job_schedule(El::LeftOrRight side, El::UpperOrLower uplo,
                           El::Orientation orientation, El::UnitOrNonUnit diag,
                           const size_t num_ranks, const size_t num_primes,
                           const std::vector<size_t> &heights,
                           const std::vector<size_t> &widths,
                           const Verbosity verbosity)
  {
    auto create_jobs_for_prime = [&](size_t prime_index, size_t split_factor) {
      ASSERT(heights.size() == widths.size());
      const auto num_blocks = heights.size();
      std::vector<Blas_Job> jobs;
      for(size_t block_index = 0; block_index < num_blocks; ++block_index)
        {
          ASSERT(split_factor != 0);
          const auto height = heights.at(block_index);
          const auto width = widths.at(block_index);

          El::IR all_I(0, height);
          El::IR all_J(0, width);
          std::vector<El::IR> Is;
          std::vector<El::IR> Js;
          if(side == El::LEFT)
            {
              // split vertically
              Is = {all_I};
              Js = split_range(all_J, split_factor);
            }
          else
            {
              // split horizontally
              Is = split_range(all_I, split_factor);
              Js = {all_J};
            }

          for(const auto I : Is)
            {
              for(const auto J : Js)
                {
                  jobs.push_back(Blas_Job::create_job(side, uplo, orientation,
                                                      diag, prime_index,
                                                      block_index, I, J));
                }
            }
        }
      return jobs;
    };

    size_t split_factor;
    return ::create_blas_job_schedule<Blas_Job>(
      num_ranks, num_primes, create_jobs_for_prime, verbosity, split_factor);
  }
}