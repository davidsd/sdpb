#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "blas_jobs/Blas_Job_Schedule.hxx"
#include "fmpz/fmpz_mul_blas_util.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <cblas.h>

// code adopted from flint mul_blas.c
namespace
{
  El::Int sum(const std::vector<El::Int> &block_heights)
  {
    return std::accumulate(block_heights.begin(), block_heights.end(), 0);
  }

  template <class T>
  size_t window_size_bytes(const Residue_Matrices_Window<T> &window)
  {
    return window.height * window.width * window.num_primes * sizeof(T);
  }
}

BigInt_Shared_Memory_Syrk_Context::BigInt_Shared_Memory_Syrk_Context(
  const El::mpi::Comm &shared_memory_comm, size_t group_index,
  mp_bitcnt_t precision, size_t max_shared_memory_bytes,
  const std::vector<El::Int> &blocks_height_per_group, int block_width,
  const std::vector<size_t> &block_index_local_to_global, bool debug,
  const std::function<Blas_Job_Schedule(size_t num_ranks, size_t num_primes,
                                        int output_width, bool debug)>
    &create_job_schedule)
    : shared_memory_comm(shared_memory_comm),
      group_index(group_index),
      num_groups(blocks_height_per_group.size()),
      comb(precision, precision, 1, sum(blocks_height_per_group)),
      output_residues_window(shared_memory_comm, comb.num_primes, block_width,
                             block_width, debug),
      block_index_local_to_global(block_index_local_to_global),
      blas_job_schedule(create_job_schedule(
        shared_memory_comm.Size(), comb.num_primes, block_width, debug))
{
  // TODO split output window too
  const size_t output_window_bytes = window_size_bytes(output_residues_window);
  ASSERT(max_shared_memory_bytes > output_window_bytes,
         DEBUG_STRING(max_shared_memory_bytes),
         DEBUG_STRING(output_window_bytes));

  const size_t max_input_window_bytes
    = max_shared_memory_bytes - output_window_bytes;
  const size_t total_block_height = sum(blocks_height_per_group);
  size_t total_input_residues_bytes
    = total_block_height * block_width * comb.num_primes * sizeof(double);
  if(total_input_residues_bytes <= max_input_window_bytes)
    {
      // no need to split
      input_grouped_block_residues_window
        = std::make_unique<Block_Residue_Matrices_Window<double>>(
          shared_memory_comm, comb.num_primes, blocks_height_per_group.size(),
          blocks_height_per_group, block_width, debug);
    }
  else
    {
      // TODO split
      RUNTIME_ERROR("Cannot allocated shared memory window for block "
                    "residues, required: ",
                    total_input_residues_bytes,
                    " bytes, available: ", max_input_window_bytes, " bytes, ",
                    DEBUG_STRING(max_shared_memory_bytes),
                    DEBUG_STRING(output_window_bytes));
    }

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
        = input_grouped_block_residues_window->residues.at(job.prime_index);
      auto &output_matrix
        = output_residues_window.residues.at(job.prime_index);

      El::Range<El::Int> all_rows(0, input_matrix.Height());
      El::Matrix<double> submatrix;

      // P_I = 0
      El::View(submatrix, input_matrix, all_rows, job.I);
      ASSERT_EQUAL(submatrix.LockedBuffer(),
                   input_matrix.LockedBuffer(0, job.I.beg));
      El::Zero(submatrix);

      // P_J = 0
      if(job.I.beg != job.J.beg)
        {
          El::View(submatrix, input_matrix, all_rows, job.J);
          ASSERT_EQUAL(submatrix.LockedBuffer(),
                       input_matrix.LockedBuffer(0, job.J.beg));
          El::Zero(submatrix);
        }

      // Q_IJ = 0
      El::View(submatrix, output_matrix, job.I, job.J);
      ASSERT_EQUAL(submatrix.LockedBuffer(),
                   output_matrix.LockedBuffer(job.I.beg, job.J.beg));
      El::Zero(submatrix);
    }
  input_grouped_block_residues_window->Fence();
  output_residues_window.Fence();
}

El::Int BigInt_Shared_Memory_Syrk_Context::input_group_height_per_prime() const
{
  return input_grouped_block_residues_window->block_residues.at(0)
    .at(group_index)
    .Height();
}
