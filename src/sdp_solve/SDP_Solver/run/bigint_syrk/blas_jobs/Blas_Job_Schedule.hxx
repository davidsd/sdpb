#pragma once

#include "Blas_Job.hxx"

// Blas_Job_Schedule distributes BLAS jobs to ranks on a node.
// Used by BigInt_Shared_Memory_Syrk_Context
struct Blas_Job_Schedule
{
  const std::vector<std::vector<Blas_Job>> jobs_by_rank;

  Blas_Job_Schedule(size_t num_ranks, const std::vector<Blas_Job> &jobs);

  // Max cost among all ranks
  [[nodiscard]] Blas_Job::Cost max_rank_cost() const;
};
