#include "create_blas_job_schedule.hxx"

#include "sdpb_util/bigint_shared_memory/blas_jobs/create_blas_job_schedule.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/split_range.hxx"

using Job_Schedule = Blas_Job_Schedule<Blas_Job>;

std::vector<Blas_Job> create_syrk_or_gemm_jobs_for_prime(
  const Blas_Job::Kind kind, const El::UpperOrLower uplo,
  const El::Int output_matrix_height, const El::Int output_matrix_width,
  const size_t prime_index, const size_t split_factor)
{
  std::vector<Blas_Job> jobs;
  const auto ranges_I = split_range({0, output_matrix_height}, split_factor);
  const auto ranges_J = split_range({0, output_matrix_width}, split_factor);
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
  return jobs;
}

// Used in unit tests
Job_Schedule create_syrk_or_gemm_job_schedule_split_remaining_primes(
  const Blas_Job::Kind kind, const El::UpperOrLower uplo,
  const size_t num_ranks, const size_t num_primes,
  const El::Int output_matrix_height, const El::Int output_matrix_width,
  const size_t split_factor)
{
  return create_blas_job_schedule_split_remaining_primes<Blas_Job>(
    num_ranks, num_primes,
    [&kind, &uplo, &output_matrix_height, &output_matrix_width](
      const size_t prime_index, const size_t split_factor) {
      return create_syrk_or_gemm_jobs_for_prime(
        kind, uplo, output_matrix_height, output_matrix_width, prime_index,
        split_factor);
    },
    split_factor);
}

// Jobs for C = A^T B modulo each prime
// If A==B, set kind=Syrk, otherwise kind=Gemm
Job_Schedule
create_blas_job_schedule(const Blas_Job::Kind kind,
                         const El::UpperOrLower uplo, const size_t num_ranks,
                         const size_t num_primes,
                         const El::Int output_matrix_height,
                         const El::Int output_matrix_width,
                         const Verbosity verbosity)
{
  auto create_jobs_for_prime = [&](const size_t prime_index,
                                   const size_t split_factor) {
    return create_syrk_or_gemm_jobs_for_prime(kind, uplo, output_matrix_height,
                                              output_matrix_width, prime_index,
                                              split_factor);
  };

  size_t split_factor;
  const auto result = create_blas_job_schedule<Blas_Job>(
    num_ranks, num_primes, create_jobs_for_prime, verbosity, split_factor);
  if(verbosity >= Verbosity::trace && El::mpi::Rank() % num_ranks == 0)
    {
      std::ostringstream ss;
      El::BuildStream(ss,
                      "\n\tChoosing the optimal split_factor=", split_factor);
      El::BuildStream(ss, "\n\tFor the first ",
                      num_primes - num_primes % num_ranks,
                      " primes p_k, (Q mod p_k) will be calculated"
                      " via single cblas_dsyrk() call.");
      El::BuildStream(ss, "\n\tFor each of the remaining (", num_primes,
                      " mod ", num_ranks, ") = ", num_primes % num_ranks,
                      " primes, (Q mod p_k) will be split into ", split_factor,
                      " diagonal blocks and ",
                      split_factor * (split_factor - 1) / 2,
                      " above-diagonal blocks.");
      El::BuildStream(ss,
                      "\n\tEach block is calculated independently"
                      " via cblas_dsyrk() or cblas_dgemm().\n"
                      "The jobs are distributed across all ",
                      num_ranks, " ranks on the node.");
      El::BuildStream(ss, "\n\t-----------------------------");
      El::Output(ss.str());
    }
  return result;
}
