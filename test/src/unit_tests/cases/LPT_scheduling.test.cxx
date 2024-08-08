#include <algorithm>
#include "sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/LPT_scheduling.hxx"
#include "unit_tests/util/util.hxx"

namespace
{
  struct Simple_Job
  {
    const size_t id;
    const size_t cost;

    Simple_Job(const size_t id, const size_t cost) : id(id), cost(cost) {}
  };

  using Scheduling = std::vector<std::vector<Simple_Job>>;

  // Assign all jobs to the first process
  Scheduling single_process_scheduling(size_t num_ranks,
                                       const std::vector<Simple_Job> &jobs)
  {
    Scheduling jobs_by_rank(num_ranks);
    for(const auto &job : jobs)
      {
        jobs_by_rank.at(0).emplace_back(job);
      }
    return jobs_by_rank;
  }

  // Naive round-robin scheduling, ignoring job costs
  Scheduling
  round_robin_scheduling(size_t num_ranks, const std::vector<Simple_Job> &jobs)
  {
    Scheduling jobs_by_rank(num_ranks);
    for(size_t i = 0; i < jobs.size(); ++i)
      {
        size_t rank = i % num_ranks;
        jobs_by_rank.at(rank).emplace_back(jobs.at(i));
      }
    return jobs_by_rank;
  }

  // maximum completion time among all ranks
  size_t max_total_cost(const Scheduling &jobs_by_rank)
  {
    size_t max_total_cost = 0;
    for(const auto &jobs : jobs_by_rank)
      {
        size_t total_cost = 0;
        for(const auto &job : jobs)
          total_cost += job.cost;
        max_total_cost = std::max(total_cost, max_total_cost);
      }
    return max_total_cost;
  }

  // Tests that LPT works and is good enough.
  // Returns (LPT_max_cost, round_robin_max_cost, single_process_max_cost)
  // where XXX_max_cost is maximal total cost of a rank.
  std::tuple<size_t, size_t, size_t>
  test_scheduling(size_t num_ranks, std::vector<size_t> costs)
  {
    size_t num_jobs = costs.size();
    size_t total_cost = std::accumulate(costs.begin(), costs.end(), (size_t)0);

    CAPTURE(num_ranks);
    CAPTURE(num_jobs);
    CAPTURE(total_cost);
    CAPTURE(costs);

    std::vector<Simple_Job> jobs;
    for(size_t i = 0; i < num_jobs; ++i)
      jobs.emplace_back(i, costs.at(i));

    std::function<size_t(const Simple_Job &)> get_cost
      = [](const auto &job) { return job.cost; };

    auto LPT_result = LPT_scheduling(num_ranks, jobs, get_cost);
    auto round_robin_result = round_robin_scheduling(num_ranks, jobs);
    auto single_process_result = single_process_scheduling(num_ranks, jobs);

    auto LPT_max_cost = max_total_cost(LPT_result);
    auto round_robin_max_cost = max_total_cost(round_robin_result);
    auto single_process_max_cost = max_total_cost(single_process_result);
    CAPTURE(LPT_max_cost);
    CAPTURE(round_robin_max_cost);
    CAPTURE(single_process_max_cost);
    REQUIRE(single_process_max_cost == total_cost);

    // Print all costs
    std::vector<std::vector<size_t>> LPT_costs_by_rank(num_ranks);
    std::vector<std::vector<size_t>> round_robin_costs_by_rank(num_ranks);
    size_t total_LPT_cost = 0;
    size_t total_round_robin_cost = 0;
    for(size_t rank = 0; rank < num_ranks; ++rank)
      {
        for(auto &job : LPT_result.at(rank))
          {
            LPT_costs_by_rank.at(rank).push_back(job.cost);
            total_LPT_cost += job.cost;
          }
        std::sort(LPT_costs_by_rank.at(rank).begin(),
                  LPT_costs_by_rank.at(rank).end());
        for(auto &job : round_robin_result.at(rank))
          {
            round_robin_costs_by_rank.at(rank).push_back(job.cost);
            total_round_robin_cost += job.cost;
          }
        std::sort(round_robin_costs_by_rank.at(rank).begin(),
                  round_robin_costs_by_rank.at(rank).end());
      }
    REQUIRE(total_LPT_cost == total_cost);
    REQUIRE(total_round_robin_cost == total_cost);
    CAPTURE(LPT_costs_by_rank);
    CAPTURE(round_robin_costs_by_rank);

    if(num_ranks == 1 || num_jobs == 1)
      {
        REQUIRE(LPT_max_cost == round_robin_max_cost);
        REQUIRE(round_robin_max_cost == single_process_max_cost);
      }
    else
      {
        // In the worst case, LPT result is worse than optimal
        // by a factor (4/3 - 1/3m), where m = num_jobs.
        // It is possible that round-robin gives a better result.
        // E.g., for two ranks and costs={6,8,5,4,7}
        // LPT gives {8,5,4},{7,6} (LPT_cost=17),
        // round-robin - {6,5,4},{8,7} (round_robin_cost=15).
        // Thus, we cannot require LPT_cost < round_robin_cost.
        // Instead, we require LPT_cost < round_robin_cost * LPT_factor.
        double LPT_factor = (4.0 / 3 - 1.0 / 3 / num_ranks);
        REQUIRE((double)LPT_max_cost
                <= (double)round_robin_max_cost * LPT_factor + 0.0001);
        REQUIRE(round_robin_max_cost < single_process_max_cost);
      }
    return std::tie(LPT_max_cost, round_robin_max_cost,
                    single_process_max_cost);
  }
}

TEST_CASE("LPT_scheduling")
{
  // this test is purely single-process, thus testing at rank=0 is sufficient
  if(El::mpi::Rank() != 0)
    return;

  INFO("Test Longest-processing-time-first (LPT) scheduling.");
  INFO("It is a greedy identical-machine scheduling algorithm, "
       "minimizing the maximum completion time.");
  INFO(
    "https://en.wikipedia.org/wiki/Longest-processing-time-first_scheduling");
  INFO("We generate jobs with different costs and check that LPT "
       "works good enough, comparing it to a simple round-robin algorithm.");

  SECTION("Deterministic bad cases")
  {
    SECTION("num_ranks=2, costs={6, 8, 5, 7, 4}")
    {
      INFO("Special case: round-robin is optimal and LPT is not");
      size_t num_ranks = 2;
      std::vector<size_t> costs = {6, 8, 5, 7, 4};

      auto [LPT_max_cost, round_robin_max_cost, single_process_max_cost]
        = test_scheduling(num_ranks, costs);

      REQUIRE(LPT_max_cost == 17);            // {8,5,4}, {6,7}
      REQUIRE(round_robin_max_cost == 15);    // {6,5,4}, {8,7}
      REQUIRE(single_process_max_cost == 30); // {6,8,5,7,4}, {}
    }

    SECTION("Worst case for LPT")
    {
      INFO("Checking the worst case for LPT");
      INFO("https://en.wikipedia.org/wiki/"
           "Longest-processing-time-first_scheduling#Worst-case_maximum_sum");
      size_t num_ranks = GENERATE(2, 4, 6, 8, 10);
      DYNAMIC_SECTION("num_ranks=" << num_ranks)
      {
        std::vector<size_t> costs;
        size_t m = num_ranks;

        // Generate costs in such order that round-robin
        // gives an optimal result for m ranks:
        // m,m,m
        // 2m-1, m+1
        // 2m-1, m+1
        // 2m-2, m+2
        // 2m-2, m+2
        // ...
        // 3m/2+1, 3m/2-1
        // 3m/2, 3m/2
        {
          costs.push_back(m);
          for(size_t cost = 2 * m - 1; cost > 3 * m / 2; --cost)
            {
              costs.push_back(cost);
              costs.push_back(cost);
            }
          costs.push_back(3 * m / 2);
          REQUIRE(costs.size() == m);

          costs.push_back(m);
          for(size_t cost = m + 1; cost < 3 * m / 2; ++cost)
            {
              costs.push_back(cost);
              costs.push_back(cost);
            }
          costs.push_back(3 * m / 2);
          costs.push_back(m);
          CAPTURE(costs);
          REQUIRE(costs.size() == 2 * m + 1);
        }

        auto [LPT_max_cost, round_robin_max_cost, single_process_max_cost]
          = test_scheduling(num_ranks, costs);

        REQUIRE(LPT_max_cost == 4 * m - 1);
        REQUIRE(round_robin_max_cost == 3 * m);
        REQUIRE(single_process_max_cost == 3 * m * m);
      }
    }
  }

  SECTION("Random")
  {
    size_t num_ranks = GENERATE(1, 2, 10, 100);
    size_t total_cost
      = GENERATE(1, 100, 1000, std::numeric_limits<size_t>::max());

    std::vector<size_t> costs = Test_Util::random_split(total_cost);
    CAPTURE(costs);
    size_t num_jobs = costs.size();

    DYNAMIC_SECTION("num_ranks=" << num_ranks)
    DYNAMIC_SECTION("num_jobs=" << num_jobs)
    {
      test_scheduling(num_ranks, costs);
    }
  }
}
