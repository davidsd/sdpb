#include <cblas.h>
#include <flint/fmpz.h>

#include <El.hpp>

#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "fmpz_mul_blas_util.hxx"
#include "fmpz_BigFloat_convert.hxx"
#include "restore_matrix_from_residues.hxx"
#include "compute_matrix_residues.hxx"

// code adopted from flint mul_blas.c
namespace
{
  El::Int sum(const std::vector<El::Int> &block_heights)
  {
    return std::accumulate(block_heights.begin(), block_heights.end(), 0);
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
}

BigInt_Shared_Memory_Syrk_Context::BigInt_Shared_Memory_Syrk_Context(
  const El::mpi::Comm &shared_memory_comm, mp_bitcnt_t precision,
  const std::vector<El::Int> &block_heights, El::Int block_width,
  const std::vector<size_t> &block_index_local_to_shmem)
    : shared_memory_comm(shared_memory_comm),
      comb(precision, precision, 1, sum(block_heights)),
      input_block_residues_window(shared_memory_comm, comb.num_primes,
                                  block_heights.size(), block_heights,
                                  block_width),
      output_residues_window(shared_memory_comm, comb.num_primes, block_width,
                             block_width),
      block_index_local_to_shmem(block_index_local_to_shmem)
{}

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
    for(size_t prime_index = 0; prime_index < comb.num_primes; ++prime_index)
      {
        // TODO this is simple round-robin, we can use something better
        // e.g. when there are more CPUs on the node than primes
        if(prime_index % El::mpi::Size(comm) == El::mpi::Rank(comm))
          {
            const auto &input_matrix
              = input_block_residues_window.residues.at(prime_index);
            auto &output_matrix
              = output_residues_window.residues.at(prime_index);
            syrk(uplo, input_matrix, output_matrix);
          }
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
