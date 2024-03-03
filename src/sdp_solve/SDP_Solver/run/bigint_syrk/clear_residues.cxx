#include "BigInt_Shared_Memory_Syrk_Context.hxx"

// Set all input and output residues to zero
void BigInt_Shared_Memory_Syrk_Context::clear_residues()
{
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