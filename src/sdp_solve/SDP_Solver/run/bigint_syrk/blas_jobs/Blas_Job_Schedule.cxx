#include "Blas_Job_Schedule.hxx"
#include "LPT_scheduling.hxx"

namespace
{
  std::vector<std::vector<Blas_Job>>
  get_jobs_by_rank(size_t num_ranks, const std::vector<Blas_Job> &jobs)
  {
    const auto get_cost
      = [](const Blas_Job &job) -> size_t { return job.cost(); };
    return LPT_scheduling<Blas_Job, size_t>(num_ranks, jobs, get_cost);
  }
}

Blas_Job_Schedule::Blas_Job_Schedule(size_t num_ranks,
                                     const std::vector<Blas_Job> &jobs)
    : jobs_by_rank(get_jobs_by_rank(num_ranks, jobs))
{}
