#include "Blas_Job_Schedule.hxx"

namespace
{
  std::vector<std::vector<Blas_Job>>
  get_jobs_by_rank(size_t num_ranks, const std::vector<Blas_Job> &jobs)
  {
    // Simple round-robin.
    // TODO: replace with proper job scheduling based on job costs
    std::vector<std::vector<Blas_Job>> jobs_by_rank(num_ranks);
    for(size_t i = 0; i < jobs.size(); ++i)
      {
        size_t rank = i % num_ranks;
        jobs_by_rank.at(rank).push_back(jobs[i]);
      }
    return jobs_by_rank;
  }
}

Blas_Job_Schedule::Blas_Job_Schedule(size_t num_ranks,
                                     const std::vector<Blas_Job> &jobs)
    : jobs_by_rank(get_jobs_by_rank(num_ranks, jobs))
{}
