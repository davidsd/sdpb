#include "catch2/catch_amalgamated.hpp"

#include "test_util/diff.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/create_blas_jobs_schedule.hxx"

#include <El.hpp>

using Test_Util::REQUIRE_Equal::diff;

TEST_CASE("create_blas_jobs_schedule")
{
  // this test is purely single-process, thus testing at rank=0 is sufficient
  if(El::mpi::Rank() != 0)
    return;

  size_t num_ranks = GENERATE(1, 2, 16, 128);
  size_t num_primes = GENERATE(15, 67, 128, 333);
  size_t Q_height = GENERATE(1, 100, 1011);

  DYNAMIC_SECTION("num_ranks=" << num_ranks << " num_primes=" << num_primes
                               << " Q_height=" << Q_height)
  {
    const auto Q_width = Q_height;
    // TODO test also gemm
    constexpr auto uplo = El::UPPER;
    auto schedule
      = create_blas_job_schedule(Blas_Job::syrk, uplo, num_ranks, num_primes,
                                 Q_height, Q_width, Verbosity::regular);
    {
      INFO("Check that schedule has correct number of ranks:");
      DIFF(schedule.jobs_by_rank.size(), num_ranks);
    }

    Blas_Job::Cost total_cost(0, 0);
    size_t num_jobs = 0;
    for(const auto &jobs : schedule.jobs_by_rank)
      for(const auto &job : jobs)
        {
          if(job.I.beg != job.J.beg)
            {
              INFO(job.I.beg);
              INFO(job.J.beg);
              REQUIRE((job.I.beg < job.J.beg) == (uplo == El::UPPER));
            }
          total_cost += job.cost;
          num_jobs++;
        }

    {
      INFO("Check that total number of matrix elements to calculate:");
      DIFF(total_cost.elements, Q_height * (Q_height + 1) / 2 * num_primes);
    }

    {
      // TODO we don't know split_factor here, so we cannot tell exactly
      //  how many blas calls are here.
      INFO("Check the number of blas_calls:");
      DIFF(total_cost.blas_calls, num_jobs);
      REQUIRE(total_cost.blas_calls >= num_primes);
    }

    {
      INFO("Check that max_cost is above the average cost:");
      REQUIRE(schedule.max_rank_cost().elements * num_ranks
              >= total_cost.elements);
    }
  }
}
