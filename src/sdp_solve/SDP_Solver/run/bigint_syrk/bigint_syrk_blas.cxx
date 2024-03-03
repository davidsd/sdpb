#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "fmpz/Fmpz_BigInt.hxx"
#include "reduce_scatter_DistMatrix.hxx"
#include "restore_matrix_from_residues.hxx"
#include "sdpb_util/assert.hxx"

#include <cblas.h>

namespace
{
  void deallocate_lower_half(El::DistMatrix<El::BigFloat> &matrix)
  {
    // Explicitly deallocate the lower half of matrix (i.e. below diagonal).
    // This significantly reduces the total amount of memory required.
    El::Matrix<El::BigFloat> &local(matrix.Matrix());
    for(int64_t row = 0; row < matrix.Height(); ++row)
      for(int64_t column = 0; column < row; ++column)
        {
          if(matrix.IsLocal(row, column))
            {
              mpf_clear(local(matrix.LocalRow(row), matrix.LocalCol(column))
                          .gmp_float.get_mpf_t());
              local(matrix.LocalRow(row), matrix.LocalCol(column))
                .gmp_float.get_mpf_t()[0]
                ._mp_d
                = nullptr;
            }
        }
  }

  // output = inputA^T * inputB
  void gemm(const El::Matrix<double> &input_A,
            const El::Matrix<double> &input_B, El::Matrix<double> &output)
  {
    // A: KxN matrix
    // A: KxM matrix
    // output = input^T * input: NxM matrix
    ASSERT_EQUAL(input_A.Height(), input_B.Height());
    ASSERT_EQUAL(output.Height(), input_A.Width());
    ASSERT_EQUAL(output.Width(), input_B.Width());

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
    const double beta = 1.0;
    double *C = output.Buffer();
    const CBLAS_INDEX ldc = output.LDim();
    // C := alpha * A^T * B + beta * C = (in our case) A^T * B + C
    cblas_dgemm(layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta,
                C, ldc);
  }

  // output = input^T * input
  void syrk(El::UpperOrLower uplo, const El::Matrix<double> &input,
            El::Matrix<double> &output)
  {
    // input: KxN matrix
    // output = input^T * input: NxN matrix
    ASSERT_EQUAL(input.Width(), output.Width());
    ASSERT_EQUAL(output.Height(), output.Width());

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo
      = uplo == El::UpperOrLowerNS::UPPER ? CblasUpper : CblasLower;
    CBLAS_TRANSPOSE Trans = CblasTrans;
    const CBLAS_INDEX N = input.Width();
    const CBLAS_INDEX K = input.Height();
    const double alpha = 1.0;
    const double *A = input.LockedBuffer();
    const CBLAS_INDEX lda = input.LDim();
    const double beta = 1.0;
    double *C = output.Buffer();
    const CBLAS_INDEX ldc = output.LDim();
    // C := alpha * A^T * A + beta * C = (in our case) A^T * A + C
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

    switch(job.kind)
      {
        case Blas_Job::syrk: {
          // take whole columns
          const auto input_matrix = El::LockedView(
            input_block_residues_window.residues.at(prime_index), El::ALL, I);

          syrk(uplo, input_matrix, output_matrix);
          break;
        }
        case Blas_Job::gemm: {
          const auto input_A = El::LockedView(
            input_block_residues_window.residues.at(prime_index), El::ALL, I);
          const auto input_B = El::LockedView(
            input_block_residues_window.residues.at(prime_index), El::ALL, J);
          gemm(input_A, input_B, output_matrix);
          break;
        }
        default: {
          El::RuntimeError("Unexpected Blas_Job::Kind=", job.kind);
        }
      }
  }

  void update_block_timings_with_syrk(
    El::Matrix<int32_t> &block_timings_ms, const Scoped_Timer &syrk_timer,
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    const std::vector<size_t> &block_index_local_to_global,
    const El::mpi::Comm &shared_memory_comm,
    const Block_Residue_Matrices_Window<double>
      &input_grouped_block_residues_window)
  {
    // Update block_timings_ms with syrk time.
    // We split total syrk time to individual blocks contributions, which are proportional to the block size.
    // To avoid double-counting, each core accounts only for the elements that it owns.
    // If block is distributed among several ranks,
    // then total block timing is a sum of contributions from all its ranks.
    constexpr auto get_size
      = [](auto &block) { return block.LocalHeight() * block.LocalWidth(); };

    // Estimate total CPU time spent on syrk for all blocks on the node,
    // assuming that time is similar for all ranks.
    const auto total_syrk_time
      = syrk_timer.elapsed_milliseconds() * shared_memory_comm.Size();
    const auto total_size = input_grouped_block_residues_window.height
                            * input_grouped_block_residues_window.width;
    // Average time spent on one block element
    const double time_per_element
      = static_cast<double>(total_syrk_time) / total_size;

    for(size_t local_block_index = 0;
        local_block_index < block_index_local_to_global.size();
        ++local_block_index)
      {
        const auto block_size
          = get_size(bigint_input_matrix_blocks.at(local_block_index));
        const double time = block_size * time_per_element;
        const auto global_block_index
          = block_index_local_to_global.at(local_block_index);
        block_timings_ms(global_block_index, 0) += std::round(time);
      }
  }
}

void BigInt_Shared_Memory_Syrk_Context::bigint_syrk_blas(
  El::UpperOrLower uplo,
  const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
  El::DistMatrix<El::BigFloat> &bigint_output, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms)
{
  // - Calculate Q = P^T P from all P blocks on the node.
  // - If we have only one node, then it's the global result (=bigint_output).
  // - If there are several nodes, we sum Q's from all nodes
  //   and write result to a global matrix bigint_output.

  Scoped_Timer timer(timers, "bigint_syrk_blas");

  // TODO: only UPPER is supported now.
  ASSERT_EQUAL(uplo, El::UPPER);

  if(El::mpi::Congruent(shared_memory_comm, bigint_output.DistComm()))
    {
      bigint_syrk_blas_shmem(uplo, bigint_input_matrix_blocks, bigint_output,
                             timers, block_timings_ms);
    }
  else
    {
      ASSERT(El::mpi::Size(shared_memory_comm)
             < El::mpi::Size(bigint_output.DistComm()));

      const El::Grid grid(shared_memory_comm);
      // TODO create once, in ctor
      El::DistMatrix<El::BigFloat> bigint_output_shmem(
        bigint_output.Height(), bigint_output.Width(), grid);
      deallocate_lower_half(bigint_output_shmem);
      bigint_syrk_blas_shmem(uplo, bigint_input_matrix_blocks,
                             bigint_output_shmem, timers, block_timings_ms);
      reduce_scatter(bigint_output, bigint_output_shmem, timers);
    }
}

// Calculate contribution to Q = P^T P from all ranks of a single node
// (i.e. single shared memory window)
void BigInt_Shared_Memory_Syrk_Context::bigint_syrk_blas_shmem(
  El::UpperOrLower uplo,
  const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
  El::DistMatrix<El::BigFloat> &bigint_output_shmem, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms)
{
  Scoped_Timer timer(timers, "shmem");
  size_t width = bigint_output_shmem.Width();

  ASSERT_EQUAL(bigint_output_shmem.Height(), bigint_output_shmem.Width());
  ASSERT_EQUAL(bigint_output_shmem.Height(), width);
  ASSERT_EQUAL(input_grouped_block_residues_window->width, width);

  ASSERT(El::mpi::Congruent(shared_memory_comm,
                            input_grouped_block_residues_window->Comm()));
  ASSERT(
    El::mpi::Congruent(shared_memory_comm, output_residues_window.Comm()));
  ASSERT(
    El::mpi::Congruent(shared_memory_comm, bigint_output_shmem.DistComm()));

  ASSERT(input_window_split_factor > 0);

  // Clear input and output windows
  clear_residues();

  // If input window is not big enough, we should fill input window
  // several times (taking different input block rows)
  // and call BLAS each time to update output window.
  for(size_t iter = 0; iter < input_window_split_factor; ++iter)
    {
      Scoped_Timer iter_timer(timers, "split_P_" + std::to_string(iter));

      // Compute block residues
      {
        // For each block group, compute residues and write them
        // to input residues window,
        // skipping first skip_rows rows
        auto skip_rows = iter * input_group_height_per_prime();
        compute_block_residues(bigint_input_matrix_blocks, skip_rows, timers,
                               block_timings_ms);
      }

      // Square each residue matrix
      {
        Scoped_Timer syrk_timer(timers, "syrk");
        {
          Scoped_Timer blas_timer(timers, "blas_jobs");
          auto shmem_rank = El::mpi::Rank(shared_memory_comm);
          for(const auto &job : blas_job_schedule.jobs_by_rank.at(shmem_rank))
            {
              do_blas_job(job, uplo, *input_grouped_block_residues_window,
                          output_residues_window);
            }
        }
        {
          Scoped_Timer fence_timer(timers, "fence");
          output_residues_window.Fence();
        }

        update_block_timings_with_syrk(
          block_timings_ms, syrk_timer, bigint_input_matrix_blocks,
          block_index_local_to_global, shared_memory_comm,
          *input_grouped_block_residues_window);
      }
    }

  // Restore bigint_output matrix from residues
  {
    Scoped_Timer restore_timer(timers, "from_residues");

    Fmpz_BigInt big_int_value;
    El::BigFloat big_float_value;
    for(int i = 0; i < bigint_output_shmem.LocalHeight(); ++i)
      for(int j = 0; j < bigint_output_shmem.LocalWidth(); ++j)
        {
          int global_i = bigint_output_shmem.GlobalRow(i);
          int global_j = bigint_output_shmem.GlobalCol(j);

          // Only half of output matrix is initialized, ignore the other one
          if(uplo == El::UpperOrLowerNS::UPPER && global_i > global_j)
            continue;
          if(uplo == El::UpperOrLowerNS::LOWER && global_i < global_j)
            continue;

          restore_matrix_from_residues(output_residues_window, global_i,
                                       global_j, comb, big_int_value);
          big_int_value.to_BigFloat(big_float_value);
          bigint_output_shmem.SetLocal(i, j, big_float_value);
        }
  }

  El::mpi::Barrier(bigint_output_shmem.DistComm());
}
