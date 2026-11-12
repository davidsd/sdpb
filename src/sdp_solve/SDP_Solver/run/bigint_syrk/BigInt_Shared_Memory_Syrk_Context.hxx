#pragma once

#include "Bigint_Syrk_Config.hxx"
#include "blas_jobs/create_blas_job_schedule.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.hxx"
#include "sdpb_util/bigint_shared_memory/Vertical_Block_Matrix_Residues_Window.hxx"
#include "sdpb_util/bigint_shared_memory/Matrix_Residues_Window.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Batch_Crt.hxx"
#include "sdpb_util/Shared_Window_Array.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <memory>
#include <optional>
#include <set>
#include <vector>

class BigInt_Shared_Memory_Syrk_Context : boost::noncopyable
{
public:
  using Job_Schedule = Blas_Job_Schedule<Blas_Job>;

  BigInt_Shared_Memory_Syrk_Context(
    const Bigint_Syrk_Config &cfg, size_t group_index,
    const std::vector<size_t> &block_index_local_to_global,
    Verbosity verbosity,
    const std::function<Job_Schedule(
      Blas_Job::Kind kind, El::UpperOrLower uplo, size_t num_ranks,
      size_t num_primes, int output_height, int output_width,
      Verbosity _verbosity)> &create_job_schedule
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
  // Faster paths of bigint_syrk_blas, all off by default; each gives the
  // same Q bit for bit. Set from Solver_Parameters (--batchResidues,
  // --residuesNoClear) and read at every call. --batchCrt and
  // --residuesStarVR, which need memory, are set by Bigint_Syrk_Config
  // instead (its batch_crt and residues_star_vr).
  struct Options
  {
    // compute_block_residues: the residues of P by BLAS (Batch_Crt::
    // residues, the same doubles as FLINT's multi-mod)
    bool batch_residues = false;
    // the residue windows are cleared only at the first call with each BLAS
    // job schedule (their first touch); later, the first pass over P
    // overwrites Q's residues (beta = 0)
    bool residues_no_clear = false;
  };
  Options options;

  void bigint_syrk_blas(El::UpperOrLower uplo,
                        const std::vector<El::DistMatrix<El::BigFloat>>
                          &bigint_input_matrix_blocks,
                        El::DistMatrix<El::BigFloat> &bigint_output,
                        Timers &timers, El::Matrix<int32_t> &block_timings_ms);

private:
  El::mpi::Comm shared_memory_comm;
  // Index of MPI group on a node
  size_t group_index;
  // Number of MPI groups on a node
  size_t num_groups;
  int total_block_height_per_node;
  Fmpz_Comb comb;
  const Verbosity verbosity;
  // All blocks from each MPI group are combined
  // into a single block in Block_Residue_Matrices_Window
  std::unique_ptr<Vertical_Block_Matrix_Residues_Window<double>>
    input_grouped_block_residues_window_A;
  std::unique_ptr<Vertical_Block_Matrix_Residues_Window<double>>
    input_grouped_block_residues_window_B;
  // How many times we should fill input window
  // to process all blocks:
  // (should be same for all ranks)
  size_t input_window_split_factor = 0;
  size_t output_window_split_factor = 0;
  std::unique_ptr<Matrix_Residues_Window<double>> output_residues_window;
  // --batchCrt (Bigint_Syrk_Config::batch_crt, and the comb's primes small
  // enough for Batch_Crt): Q_n restored once per node, column by column with
  // Batch_Crt, into a shared window of serialized BigFloats, allocated with
  // the other windows, that restore_and_reduce then copies from
  bool batch_crt_restore = false;
  // --residuesStarVR (Bigint_Syrk_Config::residues_star_vr): blocks on grids
  // with several process rows are copied to [STAR,VR] before their residues
  // are computed, so each column is written by one rank
  bool residues_star_vr = false;
  std::unique_ptr<Shared_Window_Array<El::byte>> restored_window;
  // Batch_Crt for the restore and for --batchResidues (made on first use)
  std::unique_ptr<Batch_Crt> batch_crt;
  // options.residues_no_clear: schedules whose windows were cleared (first
  // touch) once. The schedules live in blas_job_schedule_cache (a std::map,
  // never erased) as long as the context, so their addresses identify them.
  std::set<const Job_Schedule *> first_touched;
  const std::vector<size_t> block_index_local_to_global;
  std::function<Job_Schedule(Blas_Job::Kind kind, El::UpperOrLower uplo,
                             size_t num_ranks, size_t num_primes,
                             int output_height, int output_width,
                             Verbosity verbosity)>
    create_blas_job_schedule_func;
  std::map<std::tuple<Blas_Job::Kind, El::UpperOrLower, El::Int, El::Int>,
           Job_Schedule>
    blas_job_schedule_cache;

  [[nodiscard]] const Job_Schedule &
  get_blas_job_schedule(Blas_Job::Kind kind, El::UpperOrLower uplo,
                        El::Int output_height, El::Int output_width);

  void clear_residues(const Job_Schedule &blas_job_schedule);
  void compute_block_residues(
    Vertical_Block_Matrix_Residues_Window<double>
      &grouped_block_residues_window,
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    El::Int skip_rows, El::Range<El::Int> col_range, Timers &timers,
    El::Matrix<int32_t> &block_timings_ms);

  void bigint_syrk_blas_shmem_submatrix(
    El::UpperOrLower uplo,
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    const El::Range<El::Int> &output_I, const El::Range<El::Int> &output_J,
    Timers &timers, El::Matrix<int32_t> &block_timings_ms);

  void
  restore_and_reduce(std::optional<El::UpperOrLower> uplo,
                     El::DistMatrix<El::BigFloat> &output, Timers &timers);

  [[nodiscard]] El::Int input_group_height_per_prime() const;
};
