#pragma once

#include <vector>
#include "Blas_Job.hxx"

// Jobs for Q = P^T P modulo each prime
inline std::vector<Blas_Job>
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

// Create BLAS jobs for matrix multiplication Q = P^T P.
// Each job calculates some submatrix of Q modulo some prime (num_primes total).
// Since Q is symmetric, jobs are created only for the upper triangle of Q.
// Jobs are distributed among num_ranks ranks, see Blas_Job_Schedule.
// Q is a NxN matrix, N=output_matrix_height
inline std::vector<Blas_Job>
create_blas_jobs(size_t num_ranks, size_t num_primes,
                 El::Int output_matrix_height)
{
  return create_blas_jobs_no_split(num_primes, output_matrix_height);
}
