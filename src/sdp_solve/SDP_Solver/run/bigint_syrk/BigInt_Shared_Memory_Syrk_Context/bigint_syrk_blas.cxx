#include "../BigInt_Shared_Memory_Syrk_Context.hxx"
#include "../fmpz/Fmpz_BigInt.hxx"
#include "reduce_scatter_DistMatrix.hxx"
#include "restore_bigint_from_residues.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/split_range.hxx"

#include <cblas.h>

namespace
{
  // uplo == UPPER: deallocate below diagonal
  // uplo == LOWER: deallocate above diagonal
  void deallocate_unused_half(El::DistMatrix<El::BigFloat> &matrix,
                              const El::UpperOrLower uplo)
  {
    // Explicitly deallocate the lower half of matrix (i.e. below diagonal).
    // This significantly reduces the total amount of memory required.

    // Sanity check in case we pass off-diagonal block Q_IJ here
    ASSERT_EQUAL(matrix.Height(), matrix.Width());

    for(El::Int iLoc = 0; iLoc < matrix.LocalHeight(); ++iLoc)
      {
        for(El::Int jLoc = 0; jLoc < matrix.LocalWidth(); ++jLoc)
          {
            const auto i = matrix.GlobalRow(iLoc);
            const auto j = matrix.GlobalCol(jLoc);
            if(uplo == El::UPPER && i <= j)
              continue;
            if(uplo == El::LOWER && j <= i)
              continue;

            auto &gmp_float = matrix.Matrix().Ref(iLoc, jLoc).gmp_float;
            mpf_clear(gmp_float.get_mpf_t());
            gmp_float.get_mpf_t()[0]._mp_d = nullptr;
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
  void syrk(const El::UpperOrLower uplo, const El::Matrix<double> &input,
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
    const Blas_Job &job, const El::UpperOrLower uplo,
    const Block_Residue_Matrices_Window<double> &input_block_residues_window_A,
    const Block_Residue_Matrices_Window<double> &input_block_residues_window_B,
    Residue_Matrices_Window<double> &output_residues_window)
  {
    const auto prime_index = job.prime_index;
    const auto I = job.I;
    const auto J = job.J;

    auto output_matrix
      = El::View(output_residues_window.residues.at(prime_index), I, J);

    switch(job.kind)
      {
        case Blas_Job::syrk: {
          // take whole columns
          const auto input_matrix = El::LockedView(
            input_block_residues_window_A.residues.at(prime_index), El::ALL,
            I);

          syrk(uplo, input_matrix, output_matrix);
          break;
        }
        case Blas_Job::gemm: {
          const auto input_A = El::LockedView(
            input_block_residues_window_A.residues.at(prime_index), El::ALL,
            I);
          const auto input_B = El::LockedView(
            input_block_residues_window_B.residues.at(prime_index), El::ALL,
            J);
          gemm(input_A, input_B, output_matrix);
          break;
        }
        default: {
          El::RuntimeError("Unexpected Blas_Job::Kind=", job.kind);
        }
      }
  }

  void do_blas_jobs(
    const El::UpperOrLower uplo, const Blas_Job::Kind kind,
    const Blas_Job_Schedule &blas_job_schedule,
    const std::unique_ptr<Block_Residue_Matrices_Window<double>>
      &input_grouped_block_residues_window_A,
    const std::unique_ptr<Block_Residue_Matrices_Window<double>>
      &input_grouped_block_residues_window_B,
    std::unique_ptr<Residue_Matrices_Window<double>> &output_residues_window,
    const El::mpi::Comm &shared_memory_comm, Timers &timers)
  {
    // Square each residue matrix
    {
      Scoped_Timer blas_timer(timers, "blas_jobs");
      const auto shmem_rank = shared_memory_comm.Rank();
      for(const auto &job : blas_job_schedule.jobs_by_rank.at(shmem_rank))
        {
          if(kind == Blas_Job::syrk && job.I.beg != job.J.beg)
            ASSERT_EQUAL(job.I.beg < job.J.beg, uplo == El::UPPER);
          do_blas_job(job, uplo, *input_grouped_block_residues_window_A,
                      kind == Blas_Job::syrk
                        ? *input_grouped_block_residues_window_A
                        : *input_grouped_block_residues_window_B,
                      *output_residues_window);
        }
    }
    {
      Scoped_Timer fence_timer(timers, "fence");
      output_residues_window->Fence();
    }
  }

  void update_block_timings_with_syrk(
    El::Matrix<int32_t> &block_timings_ms, const Scoped_Timer &syrk_timer,
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    const std::vector<size_t> &block_index_local_to_global,
    const El::mpi::Comm &shared_memory_comm,
    const Block_Residue_Matrices_Window<double>
      &input_grouped_block_residues_window,
    const size_t total_block_height)
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
    // Total number of block elements on the node
    // NB: total_block_height > input_grouped_block_residues_window.height,
    // if split_factor > 1.
    const auto total_size
      = total_block_height * input_grouped_block_residues_window.width;
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

  ASSERT_EQUAL(bigint_output.Height(), bigint_output.Width());
  ASSERT(output_residues_window->width * output_window_split_factor
         >= bigint_output.Height());

  auto output_ranges
    = split_range({0, bigint_output.Height()}, output_window_split_factor);

  for(size_t i = 0; i < output_window_split_factor; ++i)
    {
      const auto &I = output_ranges.at(i);
      for(size_t j = i; j < output_window_split_factor; ++j)
        {
          Scoped_Timer ij_timer(timers, El::BuildString("Q_", i, "_", j));
          const auto &J = output_ranges.at(j);
          auto bigint_output_submatrix = bigint_output(I, J);
          if(El::mpi::Congruent(shared_memory_comm, bigint_output.DistComm()))
            {
              // Single-node case, no need to reduce-scatter
              bigint_syrk_blas_shmem_submatrix(
                uplo, bigint_input_matrix_blocks, bigint_output_submatrix, I,
                J, timers, block_timings_ms);
            }
          else
            {
              ASSERT(shared_memory_comm.Size()
                     < bigint_output.DistComm().Size());

              const El::Grid grid(shared_memory_comm);
              El::DistMatrix<El::BigFloat> bigint_output_shmem_submatrix(
                I.end - I.beg, J.end - J.beg, grid);
              // For diagonal blocks Q_II, we need only upper/lower half
              if(i == j)
                deallocate_unused_half(bigint_output_shmem_submatrix, uplo);

              bigint_syrk_blas_shmem_submatrix(
                uplo, bigint_input_matrix_blocks,
                bigint_output_shmem_submatrix, I, J, timers, block_timings_ms);

              // Sum the results for all nodes.
              // For diagonal Q blocks (i==j), we need only upper/lower half (for uplo==UPPER/LOWER, respectively)
              // For off-diagonal blocks, we have to synchronize all elements.
              reduce_scatter(bigint_output_submatrix, bigint_output_shmem_submatrix,
                             timers,
                             i == j ? std::make_optional(uplo) : std::nullopt);
            }
        }
    }
}

// Calculate contribution to Q_IJ = P_I^T P_J from all ranks of a single node
// (i.e. single shared memory window).
// bigint_output_shmem_submatrix is Q_IJ on the node
void BigInt_Shared_Memory_Syrk_Context::bigint_syrk_blas_shmem_submatrix(
  El::UpperOrLower uplo,
  const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
  El::DistMatrix<El::BigFloat> &bigint_output_shmem_submatrix,
  const El::Range<El::Int> &output_I, const El::Range<El::Int> &output_J,
  Timers &timers, El::Matrix<int32_t> &block_timings_ms)
{
  Scoped_Timer timer(timers, "shmem");

  ASSERT(El::mpi::Congruent(shared_memory_comm,
                            input_grouped_block_residues_window_A->Comm()));
  ASSERT(
    El::mpi::Congruent(shared_memory_comm, output_residues_window->Comm()));
  ASSERT(El::mpi::Congruent(shared_memory_comm,
                            bigint_output_shmem_submatrix.DistComm()));

  ASSERT(input_window_split_factor > 0);

  const auto output_height = output_I.end - output_I.beg;
  const auto output_width = output_J.end - output_J.beg;

  ASSERT_EQUAL(output_height, bigint_output_shmem_submatrix.Height());
  ASSERT_EQUAL(output_width, bigint_output_shmem_submatrix.Width());

  ASSERT(input_grouped_block_residues_window_A->width >= output_height,
         DEBUG_STRING(input_grouped_block_residues_window_A->width),
         DEBUG_STRING(output_height));
  ASSERT(input_grouped_block_residues_window_A->width >= output_width,
         DEBUG_STRING(input_grouped_block_residues_window_A->width),
         DEBUG_STRING(output_width));

  // syrk for diagonal blocks of Q (Q_II = P_I^T P_I),
  // gemm for off-diagonal (Q_IJ = P_I^T P_J)
  const auto kind = output_I == output_J ? Blas_Job::syrk : Blas_Job::gemm;
  const auto blas_job_schedule
    = get_blas_job_schedule(kind, uplo, output_height, output_width);

  // Clear input and output windows
  clear_residues(*blas_job_schedule);

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
        const auto skip_rows = iter * input_group_height_per_prime();
        compute_block_residues(*input_grouped_block_residues_window_A,
                               bigint_input_matrix_blocks, skip_rows, output_I,
                               timers, block_timings_ms);
        if(output_I.beg != output_J.beg)
          {
            ASSERT(output_I.end != output_J.end);
            compute_block_residues(*input_grouped_block_residues_window_B,
                                   bigint_input_matrix_blocks, skip_rows,
                                   output_J, timers, block_timings_ms);
          }
        // TODO if input_split_factor == 1, we can reuse
        // same P_I for all P_J
      }

      // Square each residue matrix
      {
        Scoped_Timer syrk_timer(timers, "syrk");
        do_blas_jobs(uplo, kind, *blas_job_schedule,
                     input_grouped_block_residues_window_A,
                     input_grouped_block_residues_window_B,
                     output_residues_window, shared_memory_comm, timers);
        update_block_timings_with_syrk(
          block_timings_ms, syrk_timer, bigint_input_matrix_blocks,
          block_index_local_to_global, shared_memory_comm,
          *input_grouped_block_residues_window_A, total_block_height_per_node);
      }
    }

  // Restore bigint_output matrix from residues
  {
    Scoped_Timer restore_timer(timers, "from_residues");

    Fmpz_BigInt big_int_value;
    El::BigFloat big_float_value;
    for(int i = 0; i < bigint_output_shmem_submatrix.LocalHeight(); ++i)
      for(int j = 0; j < bigint_output_shmem_submatrix.LocalWidth(); ++j)
        {
          const int global_i = bigint_output_shmem_submatrix.GlobalRow(i);
          const int global_j = bigint_output_shmem_submatrix.GlobalCol(j);

          if(kind == Blas_Job::syrk)
            {
              // Only half of output matrix is initialized, ignore the other one
              if(uplo == El::UpperOrLowerNS::UPPER && global_j < global_i)
                continue;
              if(uplo == El::UpperOrLowerNS::LOWER && global_i < global_j)
                continue;
            }

          restore_bigint_from_residues(*output_residues_window, global_i,
                                       global_j, comb, big_int_value);
          big_int_value.to_BigFloat(big_float_value);
          bigint_output_shmem_submatrix.SetLocal(i, j, big_float_value);
        }
  }

  El::mpi::Barrier(bigint_output_shmem_submatrix.DistComm());
}
