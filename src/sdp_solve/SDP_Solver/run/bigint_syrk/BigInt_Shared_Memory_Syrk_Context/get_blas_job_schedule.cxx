#include "../BigInt_Shared_Memory_Syrk_Context.hxx"

std::shared_ptr<Blas_Job_Schedule>
BigInt_Shared_Memory_Syrk_Context::get_blas_job_schedule(Blas_Job::Kind kind,
                                                         El::UpperOrLower uplo,
                                                         El::Int output_height,
                                                         El::Int output_width)
{
  const auto key = std::tie(kind, uplo, output_height, output_width);
  if(blas_job_schedule_cache.find(key) != blas_job_schedule_cache.end())
    return blas_job_schedule_cache.at(key);

  auto [it, res] = blas_job_schedule_cache.emplace(
    key, std::make_shared<Blas_Job_Schedule>(create_blas_job_schedule_func(
           kind, uplo, shared_memory_comm.Size(), comb.num_primes,
           output_height, output_width, verbosity)));
  ASSERT(res);
  return it->second;
}