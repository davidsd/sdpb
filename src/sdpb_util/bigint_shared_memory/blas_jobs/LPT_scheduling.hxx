#pragma once

#include <vector>
#include <functional>
#include <numeric>
#include <queue>
#include <algorithm>

// Input: number of ranks, costs for each job.
// Output: vector of vectors of job IDs for each rank.
// See https://en.wikipedia.org/wiki/Longest-processing-time-first_scheduling
// This is a greedy algorithm that always takes the heaviest job
// and assigns it to the most free process.
// The goal is to minimize the largest cost sum,
// i.e., the time when the system finishes all jobs.
// In the worst case, the largest cost sum is
// (4/3 - 1/3/num_ranks) times the optimal largest sum.
template <class TJob, class TCost>
std::vector<std::vector<TJob>>
LPT_scheduling(size_t num_ranks, const std::vector<TJob> &jobs,
               const std::function<TCost(const TJob &)> &get_cost)
{
  // output
  std::vector<std::vector<TJob>> jobs_by_rank(num_ranks);

  // Job_IDs: 0,1,2,...
  std::vector<size_t> job_ids(jobs.size());
  std::iota(job_ids.begin(), job_ids.end(), 0);

  // The heaviest jobs go first
  std::sort(job_ids.begin(), job_ids.end(),
            [&jobs, &get_cost](size_t a, size_t b) {
              auto a_cost = get_cost(jobs[a]);
              auto b_cost = get_cost(jobs[b]);
              return std::tie(a_cost, a) > std::tie(b_cost, b);
            });

  // Create min-priority queue,
  // sorted by TQueueItem = (total cost, number of jobs, rank)
  using TQueueItem = std::tuple<TCost, size_t, size_t>;
  std::priority_queue<TQueueItem, std::vector<TQueueItem>, std::greater<>>
    costs_by_rank;
  for(size_t rank = 0; rank < num_ranks; ++rank)
    {
      TCost cost{};
      size_t num_jobs = 0;
      costs_by_rank.emplace(cost, num_jobs, rank);
    }

  // Put each job into the rank that has the minimal load,
  // starting from the longest jobs
  for(auto job_id : job_ids)
    {
      auto [cost, num_jobs, rank] = costs_by_rank.top();

      // Update current rank cost
      costs_by_rank.pop();
      cost += get_cost(jobs.at(job_id));
      num_jobs++;
      costs_by_rank.emplace(cost, num_jobs, rank);

      jobs_by_rank.at(rank).push_back(jobs.at(job_id));
    }

  return jobs_by_rank;
}
