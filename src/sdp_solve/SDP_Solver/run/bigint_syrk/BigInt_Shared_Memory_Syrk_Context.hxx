#pragma once

#include "blas_jobs/create_blas_jobs_schedule.hxx"
#include "blas_jobs/Blas_Job_Schedule.hxx"
#include "fmpz/Fmpz_Comb.hxx"
#include "Block_Residue_Matrices_Window.hxx"
#include "Residue_Matrices_Window.hxx"
#include "sdpb_util/Timers/Timers.hxx"

// TODO rename
struct BigInt_Shared_Memory_Syrk_Context : boost::noncopyable
{
  El::mpi::Comm shared_memory_comm;
  // Index of MPI group on a node
  size_t group_index;
  // Sizes of MPI groups on a node
  const std::vector<int> &group_comm_sizes;
  // Number of MPI groups on a node
  size_t num_groups;
  int total_block_height_per_node;
  Fmpz_Comb comb;
  // All blocks from each MPI group are combined
  // into a single block in Block_Residue_Matrices_Window
  std::unique_ptr<Block_Residue_Matrices_Window<double>>
    input_grouped_block_residues_window;
  // How many times we should fill input window
  // to process all blocks:
  // (should be same for all ranks)
  size_t input_window_split_factor = 0;
  Residue_Matrices_Window<double> output_residues_window;
  const std::vector<size_t> block_index_local_to_global;
  const Blas_Job_Schedule blas_job_schedule;

  void clear_residues();
  BigInt_Shared_Memory_Syrk_Context(
    const El::mpi::Comm &shared_memory_comm, size_t group_index,
    const std::vector<int> &group_comm_sizes, mp_bitcnt_t precision,
    size_t max_shared_memory_bytes,
    const std::vector<El::Int> &blocks_height_per_group, int block_width,
    const std::vector<size_t> &block_index_local_to_global, bool debug,
    const std::function<Blas_Job_Schedule(size_t num_ranks, size_t num_primes,
                                          int output_width, bool debug)>
      &create_job_schedule
    = create_blas_job_schedule);

  // Calculate Q := P^T P
  //
  // P and Q are distributed among shared memory communicator, context.shared_memory_comm
  // (in practice - among all processes on a single machine).
  //
  // bigint_input_matrix_blocks - horizontal bands of P matrix
  // bigint_output - Q matrix
  //
  // Each process has some blocks of P, and their indices in communicator are stored in
  // block_indices_per_shared_memory_comm (vector of the same size as bigint_input_matrix_blocks)
  //
  // Both P and Q elements should be (big) integers;
  // We calculate residues of P modulo set of primes,
  // then multiply residue matrices via BLAS,
  // and restore Q from residues using Chinese Remainder Theorem
  //
  // If you want to square arbitrary BigFloat matrix P,
  // then use Matrix_Normalizer before and after calling this bigint_syrk_blas()
  void bigint_syrk_blas(El::UpperOrLower uplo,
                        const std::vector<El::DistMatrix<El::BigFloat>>
                          &bigint_input_matrix_blocks,
                        El::DistMatrix<El::BigFloat> &bigint_output,
                        Timers &timers, El::Matrix<int32_t> &block_timings_ms);

private:
  void compute_block_residues(
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    El::Int skip_rows, Timers &timers, El::Matrix<int32_t> &block_timings_ms);

  void bigint_syrk_blas_shmem(
    El::UpperOrLower uplo,
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    El::DistMatrix<El::BigFloat> &bigint_output_shmem, Timers &timers,
    El::Matrix<int32_t> &block_timings_ms);

  [[nodiscard]] El::Int input_group_height_per_prime() const;
};
