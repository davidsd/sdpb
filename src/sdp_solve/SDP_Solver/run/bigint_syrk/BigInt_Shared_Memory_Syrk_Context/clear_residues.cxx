#include "../BigInt_Shared_Memory_Syrk_Context.hxx"

// Set all input and output residues to zero
void BigInt_Shared_Memory_Syrk_Context::clear_residues(
  const Blas_Job_Schedule &blas_job_schedule)
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
  const auto rank = shared_memory_comm.Rank();
  for(const auto &job : blas_job_schedule.jobs_by_rank.at(rank))
    {
      auto &input_window_A = *input_grouped_block_residues_window_A;
      auto &input_window_B = input_grouped_block_residues_window_B == nullptr
                               ? *input_grouped_block_residues_window_A
                               : *input_grouped_block_residues_window_B;
      auto &input_matrix_A = input_window_A.residues.at(job.prime_index);
      auto &input_matrix_B = input_window_B.residues.at(job.prime_index);
      auto &output_matrix
        = output_residues_window->residues.at(job.prime_index);

      const El::Range<El::Int> all_rows(0, input_matrix_A.Height());
      El::Matrix<double> submatrix;

      // P_I = 0
      El::View(submatrix, input_matrix_A, all_rows, job.I);
      ASSERT_EQUAL(submatrix.LockedBuffer(),
                   input_matrix_A.LockedBuffer(0, job.I.beg));
      El::Zero(submatrix);

      // P_J = 0
      if(job.I.beg != job.J.beg
         || input_matrix_A.Buffer() != input_matrix_B.Buffer())
        {
          El::View(submatrix, input_matrix_B, all_rows, job.J);
          ASSERT_EQUAL(submatrix.LockedBuffer(),
                       input_matrix_B.LockedBuffer(0, job.J.beg));
          El::Zero(submatrix);
        }

      // Q_IJ = 0
      El::View(submatrix, output_matrix, job.I, job.J);
      ASSERT_EQUAL(submatrix.LockedBuffer(),
                   output_matrix.LockedBuffer(job.I.beg, job.J.beg));
      El::Zero(submatrix);
    }

  output_residues_window->Fence();
  input_grouped_block_residues_window_A->Fence();
  if(input_grouped_block_residues_window_B != nullptr)
    input_grouped_block_residues_window_B->Fence();
}