#pragma once

#include "Blas_Job_Cost.hxx"
#include "Blas_Job_Schedule.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <El.hpp>

#include <vector>

// Minimal split that is guaranteed to speed up computations,
// even if BLAS call overhead is huge.
// This means that we split only while each rank has only one job.
template <class TJob>
size_t minimal_split_factor(
  const size_t num_ranks, const size_t num_primes,
  const std::function<std::vector<TJob>(
    size_t prime_index, size_t split_factor)> &create_jobs_for_prime)
{
  size_t prev_num_jobs = 0;
  for(size_t split_factor = 1; split_factor <= num_ranks; ++split_factor)
    {
      const auto num_jobs
        = num_primes * create_jobs_for_prime(0, split_factor).size();
      if(num_jobs >= num_ranks)
        return split_factor;
      // Split factor is too large, no more splitting is happening
      if(num_jobs == prev_num_jobs)
        return split_factor - 1;
      prev_num_jobs = num_jobs;
    }
  // TODO shall we ever get to this point?
  return num_ranks;
}

template <class TJob>
Blas_Job_Schedule<TJob> create_blas_job_schedule_split_remaining_primes(
  const size_t num_ranks, const size_t num_primes,
  const std::function<std::vector<TJob>(
    size_t prime_index, size_t split_factor)> &create_jobs_for_prime,
  const size_t split_factor)
{
  std::vector<TJob> jobs;
  auto add_jobs
    = [&](const size_t prime_index, const size_t curr_split_factor) {
        auto curr_jobs = create_jobs_for_prime(prime_index, curr_split_factor);
        jobs.insert(jobs.end(), curr_jobs.begin(), curr_jobs.end());
      };
  // Distribute uniformly as many primes as we can, without splitting Q
  const auto num_primes_no_split = num_primes - num_primes % num_ranks;
  for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
    {
      if(prime_index < num_primes_no_split)
        add_jobs(prime_index, 1);
      else
        add_jobs(prime_index, split_factor);
    }
  return Blas_Job_Schedule(num_ranks, jobs);
}

template <class TJob>
Blas_Job_Schedule<TJob> create_blas_job_schedule(
  const size_t num_ranks, const size_t num_primes,
  const std::function<std::vector<TJob>(
    size_t prime_index, size_t split_factor)> &create_jobs_for_prime,
  const Verbosity &verbosity)
{
  const bool do_print
    = verbosity == Verbosity::trace && El::mpi::Rank() % num_ranks == 0;
  std::ostringstream ss;

  const auto min_split_factor
    = minimal_split_factor(num_ranks, num_primes, create_jobs_for_prime);

  if(do_print)
    {
      El::BuildStream(
        ss, "-----------------------------", "\nCreating BLAS job schedule...",
        "\n\tRanks per node: ", num_ranks, "\n\tNumber of primes: ");
    }

  // In most cases, checking the first five split factors is good enough,
  // it usually gives almost uniform (up to several percents) load balancing.
  // Also, we don't split too much because of overhead.
  size_t best_split_factor = min_split_factor;
  Blas_Job_Cost min_cost
    = {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()};
  for(auto split_factor = min_split_factor;
      split_factor < min_split_factor + 5; ++split_factor)
    {
      auto schedule = create_blas_job_schedule_split_remaining_primes(
        num_ranks, num_primes, create_jobs_for_prime, split_factor);
      auto cost = schedule.max_rank_cost();

      if(split_factor == min_split_factor || cost < min_cost)
        {
          min_cost = cost;
          best_split_factor = split_factor;
        }
      if(do_print)
        {
          auto [max_rank_cost, num_jobs] = cost;
          El::BuildStream(ss, "\n\t\tChecking split_factor=", split_factor,
                          ". For the heaviest rank: rank_cost/avg_cost=",
                          1.0 * max_rank_cost * num_ranks
                            / schedule.total_cost.cost,
                          " num_jobs=", num_jobs);
        }
    }
  if(do_print)
    {
      El::BuildStream(
        ss, "\n\tChoosing the optimal split_factor=", best_split_factor);
      El::BuildStream(ss, "\n\tFor the first ",
                      num_primes - num_primes % num_ranks,
                      " primes p_k, (Q mod p_k) will be calculated"
                      " via single cblas_dsyrk() call.");
      El::BuildStream(ss, "\n\tFor each of the remaining (", num_primes,
                      " mod ", num_ranks, ") = ", num_primes % num_ranks,
                      " primes, (Q mod p_k) will be split into ",
                      best_split_factor, " diagonal blocks and ",
                      best_split_factor * (best_split_factor - 1) / 2,
                      " above-diagonal blocks.");
      El::BuildStream(ss,
                      "\n\tEach block is calculated independently"
                      " via cblas_dsyrk() or cblas_dgemm().\n"
                      "The jobs are distributed across all ",
                      num_ranks, " ranks on the node.");
      El::BuildStream(ss, "\n\t-----------------------------");
      El::Output(ss.str());
    }
  // TODO here we are calculating jobs and schedule once again,
  // we should reuse what we did inside the for loop.
  // But it is fast, not a bottleneck.
  return create_blas_job_schedule_split_remaining_primes(
    num_ranks, num_primes, create_jobs_for_prime, best_split_factor);
}
