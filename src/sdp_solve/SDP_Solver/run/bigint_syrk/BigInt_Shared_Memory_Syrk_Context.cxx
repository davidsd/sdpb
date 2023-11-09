#include <cassert>
#include <cblas.h>

#include <El.hpp>

#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "fmpz_mul_blas_util.hxx"
#include "blas_jobs/create_blas_jobs_schedule.hxx"

// code adopted from flint mul_blas.c
namespace
{
  El::Int sum(const std::vector<El::Int> &block_heights)
  {
    return std::accumulate(block_heights.begin(), block_heights.end(), 0);
  }
}

BigInt_Shared_Memory_Syrk_Context::BigInt_Shared_Memory_Syrk_Context(
  const El::mpi::Comm &shared_memory_comm, mp_bitcnt_t precision,
  const std::vector<El::Int> &block_heights, El::Int block_width,
  const std::vector<size_t> &block_index_local_to_shmem,
  const std::vector<size_t> &block_index_local_to_global, bool debug,
  const std::function<Blas_Job_Schedule(size_t num_ranks, size_t num_primes,
                                        int output_width, bool debug)>
    &create_job_schedule)
    : shared_memory_comm(shared_memory_comm),
      comb(precision, precision, 1, sum(block_heights)),
      input_block_residues_window(shared_memory_comm, comb.num_primes,
                                  block_heights.size(), block_heights,
                                  block_width, debug),
      output_residues_window(shared_memory_comm, comb.num_primes, block_width,
                             block_width, debug),
      block_index_local_to_shmem(block_index_local_to_shmem),
      block_index_local_to_global(block_index_local_to_global),
      blas_job_schedule(create_job_schedule(El::mpi::Size(shared_memory_comm),
                                            comb.num_primes, block_width,
                                            debug))
{
  assert(block_index_local_to_shmem.size()
         == block_index_local_to_global.size());

  // Disable BLAS threading explicitly, each rank should work single-threaded
  openblas_set_num_threads(1);

  // "First touch": each rank is accessing (and setting to zero) the parts of shared memory window
  // that it will use for BLAS jobs.
  // This should improve BLAS performance due to better memory pinning across the cores
  // and thus faster memory access.
  // See e.g.:
  // A Performance Evaluation of MPI Shared Memory Programming, DAVID KARLBOM Master's Thesis (2016)
  // https://www.diva-portal.org/smash/get/diva2:938293/FULLTEXT01.pdf
  // NB:
  // 1) This is not optimal for calculating residues.
  // 2) Memory is allocated in pages (=4096B), and submatrix size is not equal to page size.
  //    This means that pinning is far from perfect.
  auto rank = El::mpi::Rank(shared_memory_comm);
  for(const auto &job : blas_job_schedule.jobs_by_rank.at(rank))
    {
      auto &input_matrix
        = input_block_residues_window.residues.at(job.prime_index);
      auto &output_matrix
        = output_residues_window.residues.at(job.prime_index);

      El::Range<El::Int> all_rows(0, input_matrix.Height());
      El::Matrix<double> submatrix;

      // P_I = 0
      El::View(submatrix, input_matrix, all_rows, job.I);
      assert(submatrix.LockedBuffer()
             == input_matrix.LockedBuffer(0, job.I.beg));
      El::Zero(submatrix);

      // P_J = 0
      if(job.I.beg != job.J.beg)
        {
          El::View(submatrix, input_matrix, all_rows, job.J);
          assert(submatrix.LockedBuffer()
                 == input_matrix.LockedBuffer(0, job.J.beg));
          El::Zero(submatrix);
        }

      // Q_IJ = 0
      El::View(submatrix, output_matrix, job.I, job.J);
      assert(submatrix.LockedBuffer()
             == output_matrix.LockedBuffer(job.I.beg, job.J.beg));
      El::Zero(submatrix);
    }
  input_block_residues_window.Fence();
  output_residues_window.Fence();
}
