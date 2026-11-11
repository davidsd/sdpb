#pragma once

#include "LPT_scheduling.hxx"

#include <algorithm>
#include <vector>

// Distribute jobs to ranks on a node.
// Used e.g. for BLAS jobs in BigInt_Shared_Memory_Syrk_Context
// TODO generalize? Now this code has some specific assumptions (sorted ranks, TJob::cost field)
template <class TJob> struct Blas_Job_Schedule
{
  using Cost = Blas_Job_Cost;

  const std::vector<std::vector<TJob>> jobs_by_rank;
  const Cost total_cost;

  Blas_Job_Schedule(const size_t num_ranks, const std::vector<TJob> &jobs)
      : jobs_by_rank(get_jobs_by_rank(num_ranks, jobs)),
        total_cost(get_total_cost(jobs))
  {}
  Blas_Job_Schedule(const Blas_Job_Schedule &other) = default;
  Blas_Job_Schedule(Blas_Job_Schedule &&other) noexcept = default;
  Blas_Job_Schedule &operator=(const Blas_Job_Schedule &other) = default;
  Blas_Job_Schedule &operator=(Blas_Job_Schedule &&other) noexcept = default;

  // Max cost among all ranks
  [[nodiscard]] Cost max_rank_cost() const
  {
    // jobs_by_rank are sorted in constructor, so the last rank always has the maximal cost.
    auto last_rank_jobs = jobs_by_rank.crbegin();
    return get_rank_cost(*last_rank_jobs);
  }

private:
  static Cost get_rank_cost(const std::vector<TJob> &jobs)
  {
    Cost cost{};
    for(const auto &job : jobs)
      cost += job.cost;

    return cost;
  }

  std::vector<std::vector<TJob>> static get_jobs_by_rank(
    size_t num_ranks, const std::vector<TJob> &jobs)
  {
    const std::function<Cost(const TJob &)> get_job_cost
      = [](const TJob &job) { return job.cost; };
    auto jobs_by_rank = LPT_scheduling(num_ranks, jobs, get_job_cost);

    // The last rank get the most of job
    std::sort(jobs_by_rank.begin(), jobs_by_rank.end(),
              [](const std::vector<TJob> &a, const std::vector<TJob> &b) {
                return get_rank_cost(a) < get_rank_cost(b);
              });

    return jobs_by_rank;
  }

  static Cost get_total_cost(const std::vector<TJob> &jobs)
  {
    Cost cost{};
    for(const auto &job : jobs)
      {
        cost += job.cost;
      }
    return cost;
  }
};
