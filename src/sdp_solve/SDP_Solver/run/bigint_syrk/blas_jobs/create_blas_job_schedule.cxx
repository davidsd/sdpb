#include "create_blas_jobs_schedule.hxx"

#include "Blas_Job_Schedule.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/split_range.hxx"

// Jobs for C = A^T B modulo each prime
// If A==B, set kind=Syrk, otherwise kind=Gemm
std::vector<Blas_Job>
create_blas_jobs_no_split(const Blas_Job::Kind kind, const size_t num_primes,
                          const El::Int output_matrix_height,
                          const El::Int output_matrix_width)
{
  if(kind == Blas_Job::syrk)
    ASSERT_EQUAL(output_matrix_height, output_matrix_width);

  std::vector<Blas_Job> jobs;
  jobs.reserve(num_primes);
  const El::Range<El::Int> all_I(0, output_matrix_height);
  const El::Range<El::Int> all_J(0, output_matrix_width);
  for(size_t prime_index = 0; prime_index < num_primes; ++prime_index)
    {
      switch(kind)
        {
        case Blas_Job::syrk:
          jobs.emplace_back(Blas_Job::create_syrk_job(prime_index, all_I));
          break;
        case Blas_Job::gemm:
          jobs.emplace_back(
            Blas_Job::create_gemm_job(prime_index, all_I, all_J));
          break;
        default: LOGIC_ERROR("Uknown BLAS job kind=", kind);
        }
    }
  return jobs;
}

// Minimal split that is guaranteed to speed up computations,
// even if BLAS call overhead is huge.
// This means that we split only while each rank has only one job.
size_t minimal_split_factor(const Blas_Job::Kind job_kind,
                            const size_t num_ranks, const size_t num_primes,
                            const El::Int output_matrix_height,
                            const El::Int output_matrix_width)
{
  if(job_kind == Blas_Job::syrk)
    ASSERT_EQUAL(output_matrix_height, output_matrix_width);

  if(num_primes == 0)
    return 1;

  size_t split_factor = 1;

  while(true)
    {
      const size_t next_split_factor = split_factor + 1;
      if(next_split_factor > output_matrix_height)
        break;
      if(next_split_factor > output_matrix_width)
        break;

      size_t num_jobs = 0;
      switch(job_kind)
        {
        case Blas_Job::syrk:
          // Total number of syrk and gemm jobs to calculate
          // upper (or lower) triangle of C := A^T A
          // when C is split into NxN, N = next_split_factor
          num_jobs
            = num_primes * next_split_factor * (next_split_factor + 1) / 2;
          break;
        case Blas_Job::gemm:
          // NxN gemm jobs for each prime for calculating C := A^T B:
          num_jobs = num_primes * next_split_factor * next_split_factor;
          break;
        default: LOGIC_ERROR("Invalid job_kind=", job_kind);
        }
      if(num_jobs > num_ranks)
        break;

      split_factor = next_split_factor;
    }

  return split_factor;
}

std::vector<Blas_Job> create_blas_jobs_split_remaining_primes(
  const Blas_Job::Kind kind, const El::UpperOrLower uplo,
  const size_t num_ranks, const size_t num_primes,
  const El::Int output_matrix_height, const El::Int output_matrix_width,
  const size_t split_factor)
{
  if(kind == Blas_Job::syrk)
    ASSERT_EQUAL(output_matrix_height, output_matrix_width);
  // TODO should we allow too large split factors?
  // These factors can occur in unit_tests;
  // In such case, split_range will return ranges of size 1.
  // ASSERT(split_factor <= std::max(output_matrix_height, output_matrix_width),
  //        DEBUG_STRING(split_factor), DEBUG_STRING(output_matrix_height),
  //        DEBUG_STRING(output_matrix_width));

  // Distribute uniformly as many primes as we can, without splitting Q
  const auto num_primes_no_split = num_primes - (num_primes % num_ranks);
  std::vector<Blas_Job> jobs = create_blas_jobs_no_split(
    kind, num_primes_no_split, output_matrix_height, output_matrix_width);

  const auto ranges_I = split_range({0, output_matrix_height}, split_factor);
  const auto ranges_J = split_range({0, output_matrix_width}, split_factor);
  for(size_t prime_index = num_primes_no_split; prime_index < num_primes;
      ++prime_index)
    {
      // NB: ranges_I.size() can be less than split_factor,
      // if split_factor > height
      for(size_t i = 0; i < ranges_I.size(); ++i)
        for(size_t j = 0; j < ranges_J.size(); ++j)
          {
            if(kind == Blas_Job::syrk)
              {
                // ignore lower half for syrk
                if(uplo == El::UPPER && j < i)
                  continue;
                if(uplo == El::LOWER && j < i)
                  continue;
              }
            const auto &I = ranges_I.at(i);
            const auto &J = ranges_J.at(j);
            if(i == j && kind == Blas_Job::syrk)
              jobs.push_back(Blas_Job::create_syrk_job(prime_index, I));
            else
              jobs.push_back(Blas_Job::create_gemm_job(prime_index, I, J));
          }
    }
  return jobs;
}

Blas_Job_Schedule create_blas_job_schedule_split_remaining_primes(
  const Blas_Job::Kind kind, const El::UpperOrLowerNS::UpperOrLower uplo,
  const size_t num_ranks, const size_t num_primes,
  const El::Int output_matrix_height, const El::Int output_matrix_width,
  const size_t split_factor)
{
  const auto jobs = create_blas_jobs_split_remaining_primes(
    kind, uplo, num_ranks, num_primes, output_matrix_height,
    output_matrix_width, split_factor);
  return Blas_Job_Schedule(num_ranks, jobs);
}

Blas_Job_Schedule
create_blas_job_schedule(const Blas_Job::Kind kind,
                         const El::UpperOrLowerNS::UpperOrLower uplo,
                         const size_t num_ranks, const size_t num_primes,
                         const El::Int output_matrix_height,
                         const El::Int output_matrix_width,
                         const Verbosity verbosity)
{
  const auto all_ranks_cost
    = kind == Blas_Job::syrk
        ? output_matrix_height * (output_matrix_height + 1) / 2 * num_primes
        : output_matrix_height * output_matrix_width * num_primes;

  const auto min_split_factor
    = minimal_split_factor(kind, num_ranks, num_primes % num_ranks,
                           output_matrix_height, output_matrix_width);

  if(verbosity >= Verbosity::trace && El::mpi::Rank() == 0)
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
        kind, uplo, num_ranks, num_primes, output_matrix_height,
        output_matrix_width, split_factor);
      auto cost = schedule.max_rank_cost();

      if(cost < min_cost)
        {
          min_cost = cost;
          best_split_factor = split_factor;
        }

      if(verbosity >= Verbosity::trace && El::mpi::Rank() == 0)
        {
          auto [max_rank_cost, num_jobs] = cost;
          El::Output("Checking split_factor=", split_factor,
                     ". For the heaviest rank: rank_cost/avg_cost=",
                     1.0 * max_rank_cost * num_ranks / all_ranks_cost,
                     " num_jobs=", num_jobs);
        }
    }

  if(verbosity >= Verbosity::trace && El::mpi::Rank() == 0)
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
    kind, uplo, num_ranks, num_primes, output_matrix_height,
    output_matrix_width, best_split_factor);
}
