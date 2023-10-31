#include <cblas.h>
#include <flint/fmpz.h>

#include <El.hpp>

#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "fmpz_mul_blas_util.hxx"
#include "fmpz_BigFloat_convert.hxx"
#include "restore_matrix_from_residues.hxx"
#include "compute_matrix_residues.hxx"
#include "blas_jobs/create_blas_jobs_schedule.hxx"

// code adopted from flint mul_blas.c
namespace
{
  El::Int sum(const std::vector<El::Int> &block_heights)
  {
    return std::accumulate(block_heights.begin(), block_heights.end(), 0);
  }

  // output = inputA^T * inputB
  void gemm(const El::Matrix<double> &input_A,
            const El::Matrix<double> &input_B, El::Matrix<double> &output)
  {
    // A: KxN matrix
    // A: KxM matrix
    // output = input^T * input: NxM matrix
    assert(input_A.Height() == input_B.Height());
    assert(output.Height() == input_A.Width());
    assert(output.Width() == input_B.Width());

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE TransA = CblasTrans;
    CBLAS_TRANSPOSE TransB = CblasNoTrans;
    const CBLAS_INDEX M = input_A.Width();
    const CBLAS_INDEX N = input_B.Width();
    const CBLAS_INDEX K = input_A.Height();
    const double alpha = 1.0;
    const double *A = input_A.LockedBuffer();
    const double *B = input_B.LockedBuffer();
    const CBLAS_INDEX lda = input_A.LDim();
    const CBLAS_INDEX ldb = input_B.LDim();
    const double beta = 0.0;
    double *C = output.Buffer();
    const CBLAS_INDEX ldc = output.LDim();
    // C := alpha * A^T * B + beta * C = (in our case) A^T * B
    cblas_dgemm(layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta,
                C, ldc);
  }

  // output = input^T * input
  void syrk(El::UpperOrLower uplo, const El::Matrix<double> &input,
            El::Matrix<double> &output)
  {
    // input: KxN matrix
    // output = input^T * input: NxN matrix
    assert(input.Width() == output.Width());
    assert(output.Height() == output.Width());

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo
      = uplo == El::UpperOrLowerNS::UPPER ? CblasUpper : CblasLower;
    CBLAS_TRANSPOSE Trans = CblasTrans;
    const CBLAS_INDEX N = input.Width();
    const CBLAS_INDEX K = input.Height();
    const double alpha = 1.0;
    const double *A = input.LockedBuffer();
    const CBLAS_INDEX lda = input.LDim();
    const double beta = 0.0;
    double *C = output.Buffer();
    const CBLAS_INDEX ldc = output.LDim();
    // C := alpha * A^T * A + beta * C = (in our case) A^T * A
    cblas_dsyrk(layout, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
  }

  // job: calculate submatrix Q_IJ = P_I^T * P_J (module some prime)
  // I and J are column ranges of P
  void do_blas_job(
    const Blas_Job &job, El::UpperOrLower uplo,
    const Block_Residue_Matrices_Window<double> &input_block_residues_window,
    Residue_Matrices_Window<double> &output_residues_window)
  {
    auto prime_index = job.prime_index;
    auto I = job.I;
    auto J = job.J;

    if(I.beg > J.beg && uplo == El::UpperOrLowerNS::UPPER)
      std::swap(I, J);
    if(I.beg < J.beg && uplo == El::UpperOrLowerNS::LOWER)
      std::swap(I, J);

    auto output_matrix
      = El::View(output_residues_window.residues.at(prime_index), I, J);

    // Diagonal blocks: call syrk
    if(I == J)
      {
        // whole columns
        const auto input_matrix = El::LockedView(
          input_block_residues_window.residues.at(prime_index), El::ALL, I);

        syrk(uplo, input_matrix, output_matrix);
      }
    // Off-diagonal: call gemm
    else
      {
        const auto input_A = El::LockedView(
          input_block_residues_window.residues.at(prime_index), El::ALL, I);
        const auto input_B = El::LockedView(
          input_block_residues_window.residues.at(prime_index), El::ALL, J);
        gemm(input_A, input_B, output_matrix);
      }
  }
}

BigInt_Shared_Memory_Syrk_Context::BigInt_Shared_Memory_Syrk_Context(
  const El::mpi::Comm &shared_memory_comm, mp_bitcnt_t precision,
  const std::vector<El::Int> &block_heights, El::Int block_width,
  const std::vector<size_t> &block_index_local_to_shmem, bool debug,
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
      blas_job_schedule(create_job_schedule(El::mpi::Size(shared_memory_comm),
                                            comb.num_primes, block_width,
                                            debug))
{
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

void BigInt_Shared_Memory_Syrk_Context::bigint_syrk_blas(
  El::UpperOrLower uplo,
  const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
  El::DistMatrix<El::BigFloat> &bigint_output, Timers &timers)
{
  Scoped_Timer timer(timers, "bigint_syrk_blas");
  size_t width = bigint_output.Width();

  // TODO replace asserts with El::RuntimeError?
  if(bigint_input_matrix_blocks.size() != block_index_local_to_shmem.size())
    {
      El::RuntimeError("bigint_input_matrix_blocks.size()=",
                       bigint_input_matrix_blocks.size(),
                       ", block_index_local_to_shmem.size()=",
                       block_index_local_to_shmem.size());
    }
  assert(bigint_input_matrix_blocks.size()
         == block_index_local_to_shmem.size());

  assert(bigint_output.Height() == bigint_output.Width());
  assert(bigint_output.Height() == width);
  assert(input_block_residues_window.width == width);

  assert(El::mpi::Congruent(shared_memory_comm,
                            input_block_residues_window.Comm()));
  assert(
    El::mpi::Congruent(shared_memory_comm, output_residues_window.Comm()));
  assert(El::mpi::Congruent(shared_memory_comm, bigint_output.DistComm()));

  // Compute residues
  {
    Scoped_Timer compute_residues_timer(timers,
                                        "bigint_syrk_blas.compute_residues");
    for(size_t i = 0; i < bigint_input_matrix_blocks.size(); ++i)
      {
        // NB: block_indices should enumerate all blocks
        // from all ranks in current node
        size_t block_index = block_index_local_to_shmem.at(i);
        const auto &block = bigint_input_matrix_blocks.at(i);
        assert(block.Width() == width);
        compute_matrix_residues(block_index, block, comb,
                                input_block_residues_window);
      }
    // wait for all ranks to fill blocks_window
    input_block_residues_window.Fence();
  }

  // Square each residue matrix
  {
    Scoped_Timer syrk_timer(timers, "bigint_syrk_blas.syrk");
    auto comm = input_block_residues_window.Comm();
    auto rank = El::mpi::Rank(comm);
    for(const auto &job : blas_job_schedule.jobs_by_rank.at(rank))
      {
        do_blas_job(job, uplo, input_block_residues_window,
                    output_residues_window);
      }
    output_residues_window.Fence();
  }

  // Restore bigint_output matrix from residues
  {
    Scoped_Timer restore_timer(timers, "bigint_syrk_blas.from_residues");

    fmpz_t big_int_value;
    fmpz_init(big_int_value);
    El::BigFloat big_float_value;
    for(int i = 0; i < bigint_output.LocalHeight(); ++i)
      for(int j = 0; j < bigint_output.LocalWidth(); ++j)
        {
          int global_i = bigint_output.GlobalRow(i);
          int global_j = bigint_output.GlobalCol(j);

          // Only half of output matrix is initialized, ignore the other one
          if(uplo == El::UpperOrLowerNS::UPPER && global_i > global_j)
            continue;
          if(uplo == El::UpperOrLowerNS::LOWER && global_i < global_j)
            continue;

          restore_matrix_from_residues(output_residues_window, global_i,
                                       global_j, comb, big_int_value);
          fmpz_t_to_BigFloat(big_int_value, big_float_value);
          bigint_output.SetLocal(i, j, big_float_value);
        }

    fmpz_clear(big_int_value); // TODO use RAII wrapper?
  }

  El::mpi::Barrier(bigint_output.DistComm());

  // TODO assert: after synchronize_Q Q(i,i)=2^2N
}
