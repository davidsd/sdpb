#include "Blas_Job_Schedule.hxx"
#include "LPT_scheduling.hxx"

namespace
{
  Blas_Job::Cost get_rank_cost(const std::vector<Blas_Job> &jobs)
  {
    Blas_Job::Cost cost{};
    for(const auto &job : jobs)
      cost += job.cost;

    return cost;
  }

  std::vector<std::vector<Blas_Job>>
  get_jobs_by_rank(size_t num_ranks, const std::vector<Blas_Job> &jobs)
  {
    const auto get_job_cost = [](const Blas_Job &job) { return job.cost; };
    auto jobs_by_rank = LPT_scheduling<Blas_Job, Blas_Job::Cost>(
      num_ranks, jobs, get_job_cost);

    // The last rank get the most of job
    std::sort(
      jobs_by_rank.begin(), jobs_by_rank.end(),
      [](const std::vector<Blas_Job> &a, const std::vector<Blas_Job> &b) {
        return get_rank_cost(a) < get_rank_cost(b);
      });

    return jobs_by_rank;
  }
}

Blas_Job_Schedule::Blas_Job_Schedule(size_t num_ranks,
                                     const std::vector<Blas_Job> &jobs)
    : jobs_by_rank(get_jobs_by_rank(num_ranks, jobs))
{}

Blas_Job::Cost Blas_Job_Schedule::max_rank_cost() const
{
  // jobs_by_rank are sorted in constructor, so the last rank always has the maximal cost.
  auto last_rank_jobs = jobs_by_rank.crbegin();
  return get_rank_cost(*last_rank_jobs);
}
