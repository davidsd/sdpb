#include <cassert>
#include "create_blas_jobs.hxx"

// Jobs for Q = P^T P modulo each prime
std::vector<Blas_Job>
create_blas_jobs_no_split(size_t num_primes, El::Int output_matrix_height)
{
  std::vector<Blas_Job> jobs;
  jobs.reserve(num_primes);
  El::Range<El::Int> all(0, output_matrix_height);
  for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
    {
      jobs.emplace_back(prime_index, all, all);
    }
  return jobs;
}

// Minimal split that is guaranteed to speed up computations,
// even if BLAS call overhead is huge.
// This means that we split only while each rank has only one job.
size_t minimal_split_factor(size_t num_ranks, size_t num_primes,
                            El::Int output_matrix_height)
{
  if(num_primes == 0)
    return 1;

  size_t split_factor = 1;

  while(true)
    {
      size_t next_split_factor = split_factor + 1;
      if(next_split_factor > (size_t)output_matrix_height)
        break;
      // Total number of syrk and gemm jobs to calculate upper triangle of Q
      // when Q is split into NxN, N = next_split_factor
      size_t num_jobs
        = num_primes * next_split_factor * (next_split_factor + 1) / 2;

      if(num_jobs > num_ranks)
        break;

      split_factor = next_split_factor;
    }

  return split_factor;
}

std::vector<Blas_Job>
create_blas_jobs_split_remaining_primes(size_t num_ranks, size_t num_primes,
                                        El::Int output_matrix_height,
                                        size_t split_factor)
{
  assert(split_factor <= (size_t)output_matrix_height);

  // Distribute uniformly as many primes as we can, without splitting Q
  auto num_primes_no_split = num_primes - (num_primes % num_ranks);
  std::vector<Blas_Job> jobs
    = create_blas_jobs_no_split(num_primes_no_split, output_matrix_height);

  // Split as even as possible,
  // e.g. for Q=10x10 and split_factor=3
  // Q is split into blocks of height/width {4,3,3}
  std::vector<El::Range<El::Int>> ranges(split_factor);
  El::Int start = 0;
  for(size_t n = 0; n < split_factor; ++n)
    {
      El::Int dim = output_matrix_height / split_factor;
      if(n < output_matrix_height % split_factor)
        dim++;

      ranges.at(n) = {start, start + dim};
      start += dim;
    }

  for(size_t prime_index = num_primes_no_split; prime_index < num_primes;
      ++prime_index)
    {
      for(size_t i = 0; i < split_factor; ++i)
        for(size_t j = i; j < split_factor; ++j)
          {
            El::Range<El::Int> I = ranges.at(i);
            El::Range<El::Int> J = ranges.at(j);
            jobs.emplace_back(prime_index, I, J);
          }
    }
  return jobs;
}

std::vector<Blas_Job> create_blas_jobs(size_t num_ranks, size_t num_primes,
                                       El::Int output_matrix_height)
{
  auto split_factor = minimal_split_factor(num_ranks, num_primes % num_ranks,
                                           output_matrix_height);
  return create_blas_jobs_split_remaining_primes(
    num_ranks, num_primes, output_matrix_height, split_factor);
}
