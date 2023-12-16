#include "create_blas_jobs_schedule.hxx"

#include "Blas_Job_Schedule.hxx"
#include "sdpb_util/assert.hxx"

// Jobs for Q = P^T P modulo each prime
std::vector<Blas_Job>
create_blas_jobs_no_split(size_t num_primes, El::Int output_matrix_height)
{
  std::vector<Blas_Job> jobs;
  jobs.reserve(num_primes);
  El::Range<El::Int> all(0, output_matrix_height);
  for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
    {
      jobs.emplace_back(Blas_Job::create_syrk_job(prime_index, all));
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
  ASSERT(split_factor <= (size_t)output_matrix_height);

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
            auto job = i == j ? Blas_Job::create_syrk_job(prime_index, I)
                              : Blas_Job::create_gemm_job(prime_index, I, J);
            jobs.emplace_back(job);
          }
    }
  return jobs;
}

Blas_Job_Schedule
create_blas_job_schedule_split_remaining_primes(size_t num_ranks,
                                                size_t num_primes,
                                                El::Int output_matrix_height,
                                                size_t split_factor)
{
  auto jobs = create_blas_jobs_split_remaining_primes(
    num_ranks, num_primes, output_matrix_height, split_factor);
  return Blas_Job_Schedule(num_ranks, jobs);
}

Blas_Job_Schedule
create_blas_job_schedule(size_t num_ranks, size_t num_primes,
                         El::Int output_matrix_height, bool debug)
{
  auto all_ranks_cost
    = output_matrix_height * (output_matrix_height + 1) / 2 * num_primes;

  auto min_split_factor = minimal_split_factor(
    num_ranks, num_primes % num_ranks, output_matrix_height);

  if(debug && El::mpi::Rank() == 0)
    {
      El::Output("-----------------------------");
      El::Output("Creating BLAS job schedule...");
      El::Output(
        "The jobs calculate contribution to Q = P^T * P from a single node."
        " Results from different nodes are accumulated later.");
      El::Output("Q.Height(): ", output_matrix_height);
      El::Output("Ranks per node: ", num_ranks);
      El::Output("Number of primes: ", num_primes);
      El::Output("Avg cost per rank: ", 1.0 * all_ranks_cost / num_ranks);
    }

  size_t best_split_factor = min_split_factor;
  Blas_Job::Cost min_cost
    = {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()};

  // In most cases, checking the first five split factors is good enough,
  // it usually gives almost uniform (up to several percents) load balancing.
  // Also we don't split too much because of overhead of.
  for(auto split_factor = min_split_factor;
      split_factor < min_split_factor + 5; ++split_factor)
    {
      if(split_factor > (size_t)output_matrix_height)
        break;

      auto schedule = create_blas_job_schedule_split_remaining_primes(
        num_ranks, num_primes, output_matrix_height, split_factor);
      auto cost = schedule.max_rank_cost();

      if(cost < min_cost)
        {
          min_cost = cost;
          best_split_factor = split_factor;
        }

      if(debug && El::mpi::Rank() == 0)
        {
          auto [max_rank_cost, num_jobs] = cost;
          El::Output("Checking split_factor=", split_factor,
                     ". For the heaviest rank: rank_cost/avg_cost=",
                     1.0 * max_rank_cost * num_ranks / all_ranks_cost,
                     " num_jobs=", num_jobs);
        }
    }

  if(debug && El::mpi::Rank() == 0)
    {
      El::Output("Choosing the optimal split_factor=", best_split_factor);
      El::Output("For the first ", num_primes - num_primes % num_ranks,
                 " primes p_k, (Q mod p_k) will be calculated"
                 " via single cblas_dsyrk() call.");
      El::Output("For each of the remaining (", num_primes, " mod ", num_ranks,
                 ") = ", num_primes % num_ranks,
                 " primes, (Q mod p_k) will be split into ", best_split_factor,
                 " diagonal blocks and ",
                 best_split_factor * (best_split_factor - 1) / 2,
                 " above-diagonal blocks.");
      El::Output("Each block is calculated independently"
                 " via cblas_dsyrk() or cblas_dgemm().\n"
                 "The jobs are distributed across all ",
                 num_ranks, " ranks on the node.");
      El::Output("-----------------------------");
    }

  // TODO here we are calculating jobs and schedule once again,
  // we should reuse what we did inside the for loop.
  // But it is fast, not a bottleneck.
  return create_blas_job_schedule_split_remaining_primes(
    num_ranks, num_primes, output_matrix_height, best_split_factor);
}
