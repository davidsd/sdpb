#include "../BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_BigInt.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/split_range.hxx"

#include <cblas.h>

#include <unordered_map>

// blas_jobs/create_blas_job_schedule.cxx
std::vector<Blas_Job> create_syrk_or_gemm_jobs_for_prime(
  Blas_Job::Kind kind, El::UpperOrLower uplo, El::Int output_matrix_height,
  El::Int output_matrix_width, size_t prime_index, size_t split_factor);

namespace
{
  // output = inputA^T * inputB
  void gemm(const El::Matrix<double> &input_A,
            const El::Matrix<double> &input_B, El::Matrix<double> &output,
            const double beta)
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
    double *C = output.Buffer();
    const CBLAS_INDEX ldc = output.LDim();
    // C := alpha * A^T * B + beta * C = (in our case) A^T * B + C
    cblas_dgemm(layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta,
                C, ldc);
  }

  // output = input^T * input
  void syrk(const El::UpperOrLower uplo, const El::Matrix<double> &input,
            El::Matrix<double> &output, const double beta)
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
    double *C = output.Buffer();
    const CBLAS_INDEX ldc = output.LDim();
    // C := alpha * A^T * A + beta * C = (in our case) A^T * A + C
    cblas_dsyrk(layout, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
  }

  // job: calculate submatrix Q_IJ = P_I^T * P_J (module some prime)
  // I and J are column ranges of P
  void do_blas_job(const Blas_Job &job, const El::UpperOrLower uplo,
                   const Vertical_Block_Matrix_Residues_Window<double>
                     &input_block_residues_window_A,
                   const Vertical_Block_Matrix_Residues_Window<double>
                     &input_block_residues_window_B,
                   Matrix_Residues_Window<double> &output_residues_window,
                   const double beta)
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

          syrk(uplo, input_matrix, output_matrix, beta);
          break;
        }
        case Blas_Job::gemm: {
          const auto input_A = El::LockedView(
            input_block_residues_window_A.residues.at(prime_index), El::ALL,
            I);
          const auto input_B = El::LockedView(
            input_block_residues_window_B.residues.at(prime_index), El::ALL,
            J);
          gemm(input_A, input_B, output_matrix, beta);
          break;
        }
        default: {
          El::RuntimeError("Unexpected Blas_Job::Kind=", job.kind);
        }
      }
  }

  using Row_Group = BigInt_Shared_Memory_Syrk_Context::Row_Group;

  // runs (global columns) clipped to [lo, hi)
  void clip_runs(const std::vector<El::Range<El::Int>> &runs, const El::Int lo,
                 const El::Int hi, std::vector<El::Range<El::Int>> &out)
  {
    out.clear();
    for(const auto &r : runs)
      {
        const auto b = std::max(r.beg, lo);
        const auto e = std::min(r.end, hi);
        if(b < e)
          out.emplace_back(b, e);
      }
  }

  // --groupedSyrk: the job's tile of Q (mod its prime), group by group over
  // the group's column runs only. offset_I and offset_J: the first global
  // column of P in the input windows A and B (= the output submatrix's
  // first row and column). The same integers as do_blas_job: the rows and
  // columns left out are exact zeros, and every partial sum is an exact
  // integer below 2^53, as the dense sum is.
  void do_grouped_blas_job(
    const Blas_Job &job, const El::UpperOrLower uplo,
    const Vertical_Block_Matrix_Residues_Window<double> &input_window_A,
    const Vertical_Block_Matrix_Residues_Window<double> &input_window_B,
    Matrix_Residues_Window<double> &output_residues_window,
    const double beta, const std::vector<Row_Group> &groups,
    const El::Int offset_I, const El::Int offset_J, const int skip_group)
  {
    const auto &A = input_window_A.residues.at(job.prime_index);
    const auto &B = input_window_B.residues.at(job.prime_index);
    auto &Q = output_residues_window.residues.at(job.prime_index);
    if(beta == 0.0)
      {
        auto tile = El::View(Q, job.I, job.J);
        El::Zero(tile);
      }
    const El::Int I_beg = offset_I + job.I.beg, I_end = offset_I + job.I.end;
    const El::Int J_beg = offset_J + job.J.beg, J_end = offset_J + job.J.end;
    std::vector<El::Range<El::Int>> runs_I, runs_J;
    for(size_t group = 0; group < groups.size(); ++group)
      {
        if(int(group) == skip_group)
          continue;
        const auto &g = groups[group];
        const El::Range<El::Int> rows(g.row_begin, g.row_end);
        clip_runs(g.col_runs, I_beg, I_end, runs_I);
        if(runs_I.empty())
          continue;
        if(job.kind == Blas_Job::syrk)
          {
            ASSERT(job.I == job.J);
            for(size_t a = 0; a < runs_I.size(); ++a)
              {
                const auto &Ra = runs_I[a];
                const El::Range<El::Int> ca(Ra.beg - offset_I,
                                            Ra.end - offset_I);
                const auto Pa = El::LockedView(A, rows, ca);
                {
                  auto Qaa = El::View(Q, ca, ca);
                  syrk(uplo, Pa, Qaa, 1.0);
                }
                for(size_t b = a + 1; b < runs_I.size(); ++b)
                  {
                    const auto &Rb = runs_I[b];
                    const El::Range<El::Int> cb(Rb.beg - offset_I,
                                                Rb.end - offset_I);
                    const auto Pb = El::LockedView(A, rows, cb);
                    // runs are sorted: all of a lies before b
                    if(uplo == El::UPPER)
                      {
                        auto Qab = El::View(Q, ca, cb);
                        gemm(Pa, Pb, Qab, 1.0);
                      }
                    else
                      {
                        auto Qba = El::View(Q, cb, ca);
                        gemm(Pb, Pa, Qba, 1.0);
                      }
                  }
              }
          }
        else
          {
            clip_runs(g.col_runs, J_beg, J_end, runs_J);
            for(const auto &Ra : runs_I)
              {
                const El::Range<El::Int> ca(Ra.beg - offset_I,
                                            Ra.end - offset_I);
                const auto Pa = El::LockedView(A, rows, ca);
                for(const auto &Rb : runs_J)
                  {
                    const El::Range<El::Int> cb(Rb.beg - offset_J,
                                                Rb.end - offset_J);
                    const auto Pb = El::LockedView(B, rows, cb);
                    auto Qab = El::View(Q, ca, cb);
                    gemm(Pa, Pb, Qab, 1.0);
                  }
              }
          }
      }
  }

  // --groupedSyrk: multiply-adds of a job (the cost for LPT scheduling)
  Blas_Job_Cost
  grouped_job_cost(const Blas_Job &job, const std::vector<Row_Group> &groups,
                   const El::Int offset_I, const El::Int offset_J)
  {
    size_t cost = 0, calls = 0;
    std::vector<El::Range<El::Int>> runs_I, runs_J;
    for(const auto &g : groups)
      {
        const size_t K = g.row_end - g.row_begin;
        clip_runs(g.col_runs, offset_I + job.I.beg, offset_I + job.I.end,
                  runs_I);
        if(runs_I.empty())
          continue;
        size_t nI = 0;
        for(const auto &r : runs_I)
          nI += r.end - r.beg;
        if(job.kind == Blas_Job::syrk)
          {
            cost += K * nI * (nI + 1) / 2;
            calls += runs_I.size() * (runs_I.size() + 1) / 2;
          }
        else
          {
            clip_runs(g.col_runs, offset_J + job.J.beg, offset_J + job.J.end,
                      runs_J);
            size_t nJ = 0;
            for(const auto &r : runs_J)
              nJ += r.end - r.beg;
            cost += K * nI * nJ;
            calls += runs_I.size() * runs_J.size();
          }
      }
    // writing the tile (zeroing, and the output traffic of each call)
    cost += size_t(job.I.end - job.I.beg) * (job.J.end - job.J.beg);
    return {cost, calls};
  }

  void do_blas_jobs(
    const El::UpperOrLower uplo, const Blas_Job::Kind kind,
    const BigInt_Shared_Memory_Syrk_Context::Job_Schedule &blas_job_schedule,
    const std::unique_ptr<Vertical_Block_Matrix_Residues_Window<double>>
      &input_grouped_block_residues_window_A,
    const std::unique_ptr<Vertical_Block_Matrix_Residues_Window<double>>
      &input_grouped_block_residues_window_B,
    const std::unique_ptr<Matrix_Residues_Window<double>>
      &output_residues_window,
    const El::mpi::Comm &shared_memory_comm, Timers &timers,
    const double beta = 1.0, const std::vector<Row_Group> *groups = nullptr,
    const El::Int offset_I = 0, const El::Int offset_J = 0,
    const int skip_group = -1)
  {
    // Square each residue matrix
    {
      Scoped_Timer blas_timer(timers, "blas_jobs");
      const auto shmem_rank = shared_memory_comm.Rank();
      for(const auto &job : blas_job_schedule.jobs_by_rank.at(shmem_rank))
        {
          if(kind == Blas_Job::syrk && job.I.beg != job.J.beg)
            ASSERT_EQUAL(job.I.beg < job.J.beg, uplo == El::UPPER);
          if(groups != nullptr)
            {
              do_grouped_blas_job(
                job, uplo, *input_grouped_block_residues_window_A,
                kind == Blas_Job::syrk
                  ? *input_grouped_block_residues_window_A
                  : *input_grouped_block_residues_window_B,
                *output_residues_window, beta, *groups, offset_I, offset_J,
                skip_group);
              continue;
            }
          do_blas_job(job, uplo, *input_grouped_block_residues_window_A,
                      kind == Blas_Job::syrk
                        ? *input_grouped_block_residues_window_A
                        : *input_grouped_block_residues_window_B,
                      *output_residues_window, beta);
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
    const Vertical_Block_Matrix_Residues_Window<double>
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

  // --groupedSyrk: P's row groups for this call (or the dense product)
  grouped_row_perm = nullptr;
  row_groups = nullptr;
  if(options.grouped)
    compute_row_groups(bigint_input_matrix_blocks, timers);
  else
    {
      pass_row_perm.clear();
      pass_row_groups.clear();
      grouped_stats = {};
    }

  const auto output_ranges
    = split_range({0, bigint_output.Height()}, output_window_split_factor);

  for(size_t i = 0; i < output_window_split_factor; ++i)
    {
      const auto &I = output_ranges.at(i);
      for(size_t j = i; j < output_window_split_factor; ++j)
        {
          Scoped_Timer ij_timer(timers, El::BuildString("Q_", i, "_", j));
          const auto &J = output_ranges.at(j);
          auto bigint_output_submatrix = bigint_output(I, J);
          // Call BLAS to calculate residues for Q(I,J)
          bigint_syrk_blas_shmem_submatrix(uplo, bigint_input_matrix_blocks, I,
                                           J, timers, block_timings_ms);
          std::optional<El::UpperOrLower> uplo_opt;
          if(I.beg == J.beg)
            uplo_opt = uplo;
          // Restore Q from residues and reduce-scatter over all nodes
          restore_and_reduce(uplo_opt, bigint_output_submatrix, timers);
        }
    }
}

// Calculate contribution to Q_IJ = P_I^T P_J from all ranks of a single node
// (i.e. single shared memory window).
// bigint_output_shmem_submatrix is Q_IJ on the node
void BigInt_Shared_Memory_Syrk_Context::bigint_syrk_blas_shmem_submatrix(
  const El::UpperOrLower uplo,
  const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
  const El::Range<El::Int> &output_I, const El::Range<El::Int> &output_J,
  Timers &timers, El::Matrix<int32_t> &block_timings_ms)
{
  Scoped_Timer timer(timers, "shmem");

  ASSERT(El::mpi::Congruent(shared_memory_comm,
                            input_grouped_block_residues_window_A->Comm()));
  ASSERT(
    El::mpi::Congruent(shared_memory_comm, output_residues_window->Comm()));

  ASSERT(input_window_split_factor > 0);

  const auto output_height = output_I.end - output_I.beg;
  const auto output_width = output_J.end - output_J.beg;

  ASSERT(input_grouped_block_residues_window_A->width >= output_height,
         DEBUG_STRING(input_grouped_block_residues_window_A->width),
         DEBUG_STRING(output_height));
  ASSERT(input_grouped_block_residues_window_A->width >= output_width,
         DEBUG_STRING(input_grouped_block_residues_window_A->width),
         DEBUG_STRING(output_width));

  // syrk for diagonal blocks of Q (Q_II = P_I^T P_I),
  // gemm for off-diagonal (Q_IJ = P_I^T P_J)
  const auto kind = output_I == output_J ? Blas_Job::syrk : Blas_Job::gemm;
  const auto &blas_job_schedule
    = get_blas_job_schedule(kind, uplo, output_height, output_width);
  // --groupedSyrk: the same tiles, balanced by their grouped cost (per pass
  // of the input window: each pass has its own groups)
  const bool grouped = !pass_row_perm.empty();
  auto make_grouped_schedule = [&](const std::vector<Row_Group> &groups) {
      Scoped_Timer schedule_timer(timers, "grouped_schedule");
      std::vector<El::Int> key{El::Int(kind), El::Int(uplo), output_I.beg,
                               output_I.end, output_J.beg, output_J.end};
      for(const auto &g : groups)
        {
          key.push_back(g.row_end - g.row_begin);
          key.push_back(g.col_runs.size());
          for(const auto &r : g.col_runs)
            {
              key.push_back(r.beg);
              key.push_back(r.end);
            }
        }
      if(const auto it = grouped_schedule_cache.find(key);
         it != grouped_schedule_cache.end())
        {
          grouped_schedule.emplace(it->second);
          return;
        }
      const size_t num_ranks = shared_memory_comm.Size();
      const size_t num_primes = comb.num_primes;
      std::optional<Job_Schedule> best;
      size_t best_split = 0;
      // split_mode 0: as SDPB's schedule (whole Q for the first primes,
      // the remaining ones split); 1: every prime split
      for(int split_mode = 0; split_mode < 2; ++split_mode)
        for(const size_t split : {1, 2, 3, 4, 5, 6, 8, 10, 12, 16})
          {
            std::vector<Blas_Job> jobs;
            const size_t no_split
              = split_mode == 0 ? num_primes - num_primes % num_ranks : 0;
            for(size_t prime = 0; prime < num_primes; ++prime)
              {
                auto prime_jobs = create_syrk_or_gemm_jobs_for_prime(
                  kind, uplo, output_height, output_width, prime,
                  prime < no_split ? 1 : split);
                for(auto &job : prime_jobs)
                  {
                    job.cost = grouped_job_cost(job, groups, output_I.beg,
                                                output_J.beg);
                    jobs.push_back(job);
                  }
              }
            Job_Schedule schedule(num_ranks, jobs);
            if(!best || schedule.max_rank_cost() < best->max_rank_cost())
              {
                best.emplace(std::move(schedule));
                best_split = split_mode == 0 ? split : 100 + split;
              }
          }
      grouped_stats.split_factor = best_split;
      grouped_schedule_cache.emplace(key, *best);
      grouped_schedule.emplace(std::move(*best));
  };

  // Clear input and output windows
  // --residuesNoClear: the clearing is
  // also each rank's first touch of its parts of the windows (page
  // placement), which matters once. Afterwards it only zeroes ~(num_primes x
  // (P rows + Q) x N) doubles per node and iteration (~70 GB at N = 7,990):
  // compute_block_residues rewrites every input entry the BLAS jobs read
  // (rows past the blocks are filled with zeros), and the first pass over P
  // now overwrites Q's residues (beta = 0) instead of adding to zeros. The
  // sums are exact integers below 2^53, so the residues are the same.
  // This needs each output tile to be written by exactly one BLAS job per
  // pass: the jobs take all rows of P (El::ALL), i.e. the inner dimension is
  // never split. A schedule that split it would need beta = 1 for all but
  // one of a tile's jobs.
  const bool no_clear = options.residues_no_clear;
  bool cleared = true;
  if(!no_clear || first_touched.count(&blas_job_schedule) == 0)
    {
      clear_residues(blas_job_schedule);
      first_touched.insert(&blas_job_schedule);
    }
  else
    cleared = false;

  // If input window is not big enough, we should fill input window
  // several times (taking different input block rows)
  // and call BLAS each time to update output window.
  for(size_t iter = 0; iter < input_window_split_factor; ++iter)
    {
      Scoped_Timer iter_timer(timers, "split_P_" + std::to_string(iter));
      if(grouped)
        {
          grouped_row_perm = &pass_row_perm.at(iter);
          row_groups = &pass_row_groups.at(iter);
          make_grouped_schedule(*row_groups);
        }

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
        do_blas_jobs(uplo, kind,
                     grouped ? *grouped_schedule : blas_job_schedule,
                     input_grouped_block_residues_window_A,
                     input_grouped_block_residues_window_B,
                     output_residues_window, shared_memory_comm, timers,
                     (!cleared && iter == 0) ? 0.0 : 1.0,
                     grouped ? row_groups : nullptr, output_I.beg,
                     output_J.beg, options.grouped_test_skip_group);
        update_block_timings_with_syrk(
          block_timings_ms, syrk_timer, bigint_input_matrix_blocks,
          block_index_local_to_global, shared_memory_comm,
          *input_grouped_block_residues_window_A, total_block_height_per_node);
      }
    }

  grouped_row_perm = nullptr;
  row_groups = nullptr;
  output_residues_window->Fence();
}

// --groupedSyrk: find the groups of rows of P on the node that have the same
// nonzero columns, and the place of each row in the grouped input window,
// for each pass of the input window. Collective over shared_memory_comm;
// every rank computes the same result from the same shared mask.
bool BigInt_Shared_Memory_Syrk_Context::compute_row_groups(
  const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
  Timers &timers)
{
  Scoped_Timer timer(timers, "grouped_pattern");
  const auto previous = grouped_stats;
  grouped_schedule_cache.clear();
  pass_row_perm.clear();
  pass_row_groups.clear();
  grouped_schedule.reset();
  grouped_stats = {};

  // rows of P on the node (mask rows) and of the window (one pass)
  const El::Int H_full = total_block_height_per_node;
  const El::Int H = input_grouped_block_residues_window_A->height;
  const El::Int N = input_width;
  const size_t W = (N + 63) / 64;
  if(H_full == 0 || H == 0 || N == 0)
    return false;
  if(!grouped_mask_window)
    grouped_mask_window = std::make_unique<Shared_Window_Array<uint64_t>>(
      shared_memory_comm, size_t(H_full) * W);
  uint64_t *mask = grouped_mask_window->data();
  {
    const size_t total = size_t(H_full) * W;
    const size_t rank = shared_memory_comm.Rank();
    const size_t num_ranks = shared_memory_comm.Size();
    std::fill(mask + total * rank / num_ranks,
              mask + total * (rank + 1) / num_ranks, uint64_t(0));
  }
  grouped_mask_window->Fence();
  {
    // Each rank sets the bits of its local nonzero entries
    Scoped_Timer set_timer(timers, "set_bits");
    El::Int row_begin = group_full_offsets.at(group_index);
    std::vector<uint64_t> local;
    for(const auto &block : bigint_input_matrix_blocks)
      {
        ASSERT_EQUAL(block.Width(), N);
        const int height_local = block.LocalHeight();
        const int width_local = block.LocalWidth();
        if(height_local > 0 && width_local > 0)
          {
            local.assign(size_t(height_local) * W, 0);
            for(int jLoc = 0; jLoc < width_local; ++jLoc)
              {
                const El::Int col = block.GlobalCol(jLoc);
                const uint64_t bit = uint64_t(1) << (col % 64);
                const size_t word = col / 64;
                for(int iLoc = 0; iLoc < height_local; ++iLoc)
                  if(mpf_sgn(block.GetLocalCRef(iLoc, jLoc)
                               .gmp_float.get_mpf_t())
                     != 0)
                    local[size_t(iLoc) * W + word] |= bit;
              }
            for(int iLoc = 0; iLoc < height_local; ++iLoc)
              {
                const El::Int row = row_begin + block.GlobalRow(iLoc);
                ASSERT(row < H_full, DEBUG_STRING(row), DEBUG_STRING(H_full));
                for(size_t w = 0; w < W; ++w)
                  {
                    const uint64_t bits = local[size_t(iLoc) * W + w];
                    if(bits != 0)
                      __atomic_fetch_or(mask + size_t(row) * W + w, bits,
                                        __ATOMIC_RELAXED);
                  }
              }
          }
        row_begin += block.Height();
      }
    ASSERT_EQUAL(row_begin, group_full_offsets.at(group_index)
                              + group_full_heights.at(group_index));
  }
  grouped_mask_window->Fence();

  const size_t num_passes = input_window_split_factor;
  std::vector<std::vector<El::Int>> perms(num_passes);
  std::vector<std::vector<Row_Group>> all_groups(num_passes);
  double grouped_madds = 0;
  size_t max_groups = 0, max_runs = 0;
  for(size_t pass = 0; pass < num_passes; ++pass)
    {
      // mask row of each window row in this pass (-1: past the group's
      // rows, zero); the pass holds rows [pass h_g, (pass + 1) h_g) of
      // each MPI group g, h_g its window height (compute_block_residues)
      std::vector<El::Int> mask_row(H, -1);
      for(size_t g = 0; g < group_row_offsets.size(); ++g)
        {
          const El::Int h_g
            = input_grouped_block_residues_window_A->block_residues.at(0)
                .at(g)
                .Height();
          for(El::Int w = 0; w < h_g; ++w)
            {
              const El::Int r = El::Int(pass) * h_g + w;
              if(r < group_full_heights.at(g))
                mask_row[group_row_offsets[g] + w]
                  = group_full_offsets.at(g) + r;
            }
        }
      // Groups in order of first appearance; all-zero rows go last
      std::vector<int> row_group(H, -1);
      std::vector<El::Int> group_first_row, group_count;
      {
        std::unordered_map<uint64_t, std::vector<int>> index;
        for(El::Int row = 0; row < H; ++row)
          {
            if(mask_row[row] < 0)
              continue;
            const uint64_t *m = mask + size_t(mask_row[row]) * W;
            uint64_t hash = 1469598103934665603ULL;
            bool any = false;
            for(size_t w = 0; w < W; ++w)
              {
                any |= m[w] != 0;
                hash = (hash ^ m[w]) * 1099511628211ULL;
                hash ^= hash >> 29;
              }
            if(!any)
              continue;
            auto &candidates = index[hash];
            int group = -1;
            for(const int g : candidates)
              if(std::equal(m, m + W,
                            mask
                              + size_t(mask_row[group_first_row[g]]) * W))
                {
                  group = g;
                  break;
                }
            if(group < 0)
              {
                group = group_first_row.size();
                candidates.push_back(group);
                group_first_row.push_back(row);
                group_count.push_back(0);
              }
            row_group[row] = group;
            ++group_count[group];
          }
      }
      const size_t num_groups = group_first_row.size();
      auto &groups = all_groups[pass];
      groups.resize(num_groups);
      size_t num_runs = 0;
      El::Int begin = 0;
      for(size_t g = 0; g < num_groups; ++g)
        {
          groups[g].row_begin = begin;
          begin += group_count[g];
          groups[g].row_end = begin;
          const uint64_t *m = mask + size_t(mask_row[group_first_row[g]]) * W;
          El::Int nnz = 0;
          for(El::Int col = 0; col < N;)
            {
              if(((m[col / 64] >> (col % 64)) & 1) == 0)
                {
                  ++col;
                  continue;
                }
              El::Int end = col + 1;
              while(end < N && ((m[end / 64] >> (end % 64)) & 1) != 0)
                ++end;
              groups[g].col_runs.emplace_back(col, end);
              nnz += end - col;
              col = end;
            }
          num_runs += groups[g].col_runs.size();
          grouped_madds
            += double(group_count[g]) * double(nnz) * double(nnz);
        }
      max_groups = std::max(max_groups, num_groups);
      max_runs = std::max(max_runs, num_runs);
      auto &perm = perms[pass];
      perm.assign(H, 0);
      std::vector<El::Int> next(num_groups);
      for(size_t g = 0; g < num_groups; ++g)
        next[g] = groups[g].row_begin;
      El::Int zero_next = begin;
      for(El::Int row = 0; row < H; ++row)
        perm[row] = row_group[row] < 0 ? zero_next++ : next[row_group[row]]++;
      ASSERT_EQUAL(zero_next, H);
    }
  const double fraction
    = grouped_madds / (double(H_full) * double(N) * double(N));

  grouped_stats.num_groups = max_groups;
  grouped_stats.num_runs = max_runs;
  grouped_stats.fraction = fraction;
  // Not worth it: many groups (many small BLAS calls), or few zeros
  const bool active = max_groups > 0 && max_groups <= 4096 && fraction < 0.95;
  if(active)
    {
      pass_row_perm = std::move(perms);
      pass_row_groups = std::move(all_groups);
    }
  grouped_stats.active = active;

  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0
     && (previous.num_groups != grouped_stats.num_groups
         || previous.num_runs != grouped_stats.num_runs
         || std::abs(previous.fraction - grouped_stats.fraction) > 1e-3))
    El::Output("groupedSyrk (node of rank 0): ", max_groups,
               " row groups and ", max_runs, " column runs at most per pass, ",
               num_passes, " pass(es), ", H_full, " rows x ", N,
               " columns; multiply-adds ", fraction, " of the dense product",
               active ? "" : "; using the dense product");
  return active;
}
