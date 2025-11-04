#include "Initialize_P_Context.hxx"

#include "bigint_trmm/blas_jobs/create_blas_job_schedule.hxx"
#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/residues.hxx"

#include <cblas.h>
#include <optional>

namespace Sdpb::Sdpa
{
  namespace
  {
    // X := L X or X := X L
    void dtrmm(const El::LeftOrRight side, const El::UpperOrLower uplo,
               const El::Orientation orientation, const El::UnitOrNonUnit diag,
               const El::Matrix<double> &L, El::Matrix<double> &X)
    {
      ASSERT_EQUAL(L.Height(), L.Width());
      if(side == El::LEFT)
        ASSERT_EQUAL(L.Width(), X.Height(), DEBUG_STRING(X.Width()));
      else
        ASSERT_EQUAL(X.Width(), L.Height(), DEBUG_STRING(X.Height()));

      CBLAS_LAYOUT layout = CblasColMajor;
      CBLAS_SIDE Side;
      switch(side)
        {
        case El::LEFT: Side = CblasLeft; break;
        case El::RIGHT: Side = CblasRight; break;
        default: LOGIC_ERROR(DEBUG_STRING(side));
        }
      CBLAS_UPLO Uplo;
      switch(uplo)
        {
        case El::LOWER: Uplo = CblasLower; break;
        case El::UPPER: Uplo = CblasUpper; break;
        default: LOGIC_ERROR(DEBUG_STRING(uplo));
        }
      CBLAS_TRANSPOSE TransA;
      switch(orientation)
        {
        case El::NORMAL: TransA = CblasNoTrans; break;
        case El::TRANSPOSE: TransA = CblasTrans; break;
        case El::ADJOINT: TransA = CblasConjTrans; break;
        default: LOGIC_ERROR(DEBUG_STRING(orientation));
        }
      CBLAS_DIAG Diag;
      switch(diag)
        {
        case El::NON_UNIT: Diag = CblasNonUnit; break;
        case El::UNIT: Diag = CblasUnit; break;
        default: LOGIC_ERROR(DEBUG_STRING(diag));
        }
      const CBLAS_INDEX M = X.Height();
      const CBLAS_INDEX N = X.Width();
      constexpr double alpha = 1.0;
      const double *A = L.LockedBuffer();
      const CBLAS_INDEX lda = L.LDim();
      double *B = X.Buffer();
      const CBLAS_INDEX ldb = X.LDim();

      cblas_dtrmm(layout, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B,
                  ldb);
    }

    void
    do_blas_job(const Trmm::Blas_Job &job,
                const Vector_Matrix_Residues_Window<double> &L_window,
                Vector_Block_Diagonal_Matrix_Residues_Window<double> &G_window)
    {
      const auto &L
        = L_window.matrices_residues.at(job.prime_index).at(job.block_index);
      const auto &G
        = G_window.matrices_residues.at(job.prime_index).at(job.block_index);
      const auto &I = job.I;
      const auto &J = job.J;
      auto G_view = G(I, J);
      {
        ASSERT_EQUAL(L.Height(), L.Width());
        // TODO: Elemental does not check validity of views.
        // Extract this to a separate function and move to sdpb_util.
        ASSERT(job.I == El::ALL || (I.beg < I.end && I.end <= G.Height()),
               "Matrix view is out of range", DEBUG_STRING(I.beg),
               DEBUG_STRING(I.end), DEBUG_STRING(G.Height()));
        ASSERT(job.J == El::ALL || (J.beg < J.end && J.end <= G.Width()),
               "Matrix view is out of range", DEBUG_STRING(J.beg),
               DEBUG_STRING(J.end), DEBUG_STRING(G.Width()));
      }
      dtrmm(job.side, job.uplo, job.orientation, job.diag, L, G_view);
    }

    // TODO deduplicate with syrk code
    // TODO optimize for the case when the whole column is owned by a single rank (see syrk code).
    void compute_column_residues_elementwise(
      const El::DistMatrix<El::BigFloat> &matrix, const int global_col,
      const size_t prime_stride, Fmpz_Comb &comb,
      El::Matrix<double> &residues_first_prime)
    {
      if(!matrix.IsLocalCol(global_col))
        return;

      // Block submatrix for a given row range and a given column
      const auto block_column_submatrix
        = matrix(El::IR(0, matrix.Height()), global_col);
      ASSERT(block_column_submatrix.Viewing());

      // Residues of blocks on a current MPI group
      // modulo the first prime are written here

      // Part of shared memory window containing
      // residues of block_column_submatrix modulo the first prime number.
      const El::IR residue_I(0, matrix.Height());
      const El::IR residue_J(global_col);
      ASSERT(residue_I.beg < residues_first_prime.Height(),
             DEBUG_STRING(residue_I.beg),
             DEBUG_STRING(residues_first_prime.Height()));
      ASSERT(residue_I.end <= residues_first_prime.Height());
      auto first_residue_column_submatrix
        = residues_first_prime(residue_I, residue_J);
      ASSERT(first_residue_column_submatrix.Viewing());

      Fmpz_BigInt bigint_value;
      // Submatrix is single-column, thus index j is always 0.
      const int j = 0;
      const int jLoc = block_column_submatrix.LocalCol(j);
      for(int iLoc = 0; iLoc < block_column_submatrix.LocalHeight(); ++iLoc)
        {
          const int i = block_column_submatrix.GlobalRow(iLoc);
          bigint_value.from_BigFloat(
            block_column_submatrix.GetLocalCRef(iLoc, jLoc));
          double *data = first_residue_column_submatrix.Buffer(i, j);
          bigint_to_residues(bigint_value, comb, data, prime_stride);
        }
    }

    void compute_column_residues_elementwise(
      const El::DistMatrix<El::BigFloat> &block, const size_t block_index,
      const El::Int global_col, Fmpz_Comb &comb,
      Vector_Matrix_Residues_Window<double> &block_residues_window)
    {
      compute_column_residues_elementwise(
        block, global_col, block_residues_window.prime_stride, comb,
        block_residues_window.matrices_residues.at(0).at(block_index));
    }
    void compute_column_residues_elementwise(
      const El::DistMatrix<El::BigFloat> &block, const size_t block_index,
      const size_t vec_index, const El::Int global_col, Fmpz_Comb &comb,
      Vector_Block_Diagonal_Matrix_Residues_Window<double> &window)
    {
      compute_column_residues_elementwise(
        block, global_col, window.prime_stride, comb,
        window.prime_block_vec_matrices.at(0).at(block_index).at(vec_index));
    }

    void restore_distmatrix_from_residues(
      const std::optional<El::UpperOrLower> &uplo,
      const El::Matrix<double> &first_residue_matrix,
      const size_t prime_stride, Fmpz_Comb &comb,
      El::DistMatrix<El::BigFloat> &output)
    {
      const int height = first_residue_matrix.Height();
      const int width = first_residue_matrix.Width();
      ASSERT_EQUAL(height, output.Height());
      ASSERT_EQUAL(width, output.Width());

      auto skip_element = [&uplo](auto row, auto column) {
        if(uplo.has_value())
          {
            if(*uplo == El::UPPER && column < row)
              return true;
            if(*uplo == El::LOWER && row < column)
              return true;
          }
        return false;
      };

      std::vector<mp_limb_t> residues_buffer_temp;
      Fmpz_BigInt bigint_temp;
      for(int i = 0; i < height; ++i)
        for(int j = 0; j < width; ++j)
          {
            if(skip_element(i, j))
              continue;
            if(!output.IsLocal(i, j))
              continue;

            const double *first_residue
              = first_residue_matrix.LockedBuffer(i, j);
            residues_to_bigint(first_residue, prime_stride, comb,
                               residues_buffer_temp, bigint_temp);
            const auto iLoc = output.LocalRow(i);
            const auto jLoc = output.LocalCol(j);
            bigint_temp.to_BigFloat(output.Matrix()(iLoc, jLoc));
          }
    }
  }

  Initialize_P_Context::Initialize_P_Context(const Initialize_P_Config &cfg)
      : cfg(cfg),
        shared_memory_comm(cfg.shared_memory_comm),
        comb(cfg.comb()),
        L_X_inv_residues(create_L_window()),
        L_Y_residues(create_L_window()),
        G_window_view(shared_memory_comm, cfg.F_window_size())
  {}
  Vector_Block_Diagonal_Matrix_Residues_Window<double> &
  Initialize_P_Context::G_residues(
    const size_t num_primal_indices,
    const El::VerticalOrHorizontal vertical_or_horizontal)
  {
    const auto key = std::tie(num_primal_indices, vertical_or_horizontal);
    const auto it = G_window_cache.find(key);
    if(it != G_window_cache.end())
      return it->second;
    // TODO: for the last iteration, num_primal_indices can be smaller than primal_dimension_step,
    // In that case, residues window for G is smaller than G_window_view.
    // Shall we just increase prime_stride and fill the whole buffer (leaving gaps inside) instead?
    // If we do so, each prime will be kept at the same area as in all previous iterations.
    // This could help with NUMA memory pinning, if we .
    return G_window_cache
      .emplace(key, Vector_Block_Diagonal_Matrix_Residues_Window<double>(
                      G_window_view, comb.num_primes, cfg.node_block_dims,
                      num_primal_indices, vertical_or_horizontal))
      .first->second;
  }
  void Initialize_P_Context::compute_residues(
    const Block_Diagonal_Matrix &bdm,
    Vector_Matrix_Residues_Window<double> &window,
    El::Matrix<int32_t> &block_timings_ms)
  {
    ASSERT_EQUAL(bdm.blocks.size(), cfg.local_block_locations().size());
    for(size_t block = 0; block < bdm.blocks.size(); ++block)
      {
        const auto &loc = cfg.local_block_locations().at(block);
        const auto node_block_index = loc.block_index_node;
        const auto global_block_index = loc.block_index_global;
        const auto &matrix = bdm.blocks.at(block);
        Timer timer;
        for(int j = 0; j < matrix.Width(); ++j)
          {
            compute_column_residues_elementwise(matrix, node_block_index, j,
                                                comb, window);
          }
        block_timings_ms(global_block_index, 0)
          += timer.elapsed_milliseconds();
      }
    window.Fence();
  }
  void Initialize_P_Context::compute_residues(
    const std::vector<Block_Diagonal_Matrix> &bdms,
    Vector_Block_Diagonal_Matrix_Residues_Window<double> &window,
    El::Matrix<int32_t> &block_timings_ms)
  {
    ASSERT(bdms.size() <= cfg.primal_dimension_step,
           "vector<Block_Diagonal_Matrix> is longer than allowed by config",
           DEBUG_STRING(bdms.size()), DEBUG_STRING(cfg.primal_dimension_step));
    const auto num_blocks = cfg.local_block_locations().size();
    for(size_t vec_index = 0; vec_index < bdms.size(); ++vec_index)
      {
        const auto &bdm = bdms.at(vec_index);
        ASSERT_EQUAL(bdm.blocks.size(), num_blocks);
        for(size_t block = 0; block < num_blocks; ++block)
          {
            const auto &loc = cfg.local_block_locations().at(block);
            const auto node_block_index = loc.block_index_node;
            const auto global_block_index = loc.block_index_global;
            const auto &matrix = bdm.blocks.at(block);
            Timer timer;
            for(int j = 0; j < matrix.Width(); ++j)
              {
                // TODO add more asserts to check dimensions, array sizes etc.
                compute_column_residues_elementwise(
                  matrix, node_block_index, vec_index, j, comb, window);
              }
            block_timings_ms(global_block_index, 0)
              += timer.elapsed_milliseconds();
          }
      }
    window.Fence();
  }

  void Initialize_P_Context::trmm(
    const El::LeftOrRight side, const El::UpperOrLower uplo,
    const El::Orientation orientation, const El::UnitOrNonUnit diag,
    Vector_Matrix_Residues_Window<double> &L_window,
    Vector_Block_Diagonal_Matrix_Residues_Window<double> &G_window,
    std::vector<Block_Diagonal_Matrix> &G, const Verbosity verbosity,
    Timers &timers, El::Matrix<int32_t> &block_timings_ms)
  {
    Scoped_Timer timer(timers, "trmm");
    std::vector<size_t> heights;
    std::vector<size_t> widths;
    for(auto &stacked_matrix : G_window.matrices_residues.at(0))
      {
        heights.push_back(stacked_matrix.Height());
        widths.push_back(stacked_matrix.Width());
      }

    // TODO cache it
    const auto blas_job_schedule = Trmm::create_blas_job_schedule(
      side, uplo, orientation, diag, shared_memory_comm.Size(),
      cfg.num_primes(), heights, widths, verbosity);

    do_blas_jobs(blas_job_schedule, L_window, G_window, timers,
                 block_timings_ms);
    restore_from_residues(G_window, G, timers, block_timings_ms);
  }
  Vector_Matrix_Residues_Window<double>
  Initialize_P_Context::create_L_window() const
  {
    const Shared_Window_Array_View<double> win_view(shared_memory_comm,
                                                    cfg.L_window_size());
    return Vector_Matrix_Residues_Window(
      win_view, comb.num_primes, cfg.node_block_dims, cfg.node_block_dims);
  }
  void Initialize_P_Context::do_blas_jobs(
    const Blas_Job_Schedule<Trmm::Blas_Job> &schedule,
    Vector_Matrix_Residues_Window<double> &L_window,
    Vector_Block_Diagonal_Matrix_Residues_Window<double> &G_window,
    Timers &timers, El::Matrix<int32_t> &block_timings_ms)
  {
    {
      Scoped_Timer blas_timer(timers, "blas_jobs");
      const auto shmem_rank = shared_memory_comm.Rank();
      for(const auto &job : schedule.jobs_by_rank.at(shmem_rank))
        {
          const auto &loc = cfg.node_block_locations().at(job.block_index);
          ASSERT_EQUAL(loc.block_index_node, job.block_index);
          const auto block_index = loc.block_index_global;
          Timer timer;
          do_blas_job(job, L_window, G_window);
          block_timings_ms(block_index, 0) += timer.elapsed_milliseconds();
        }
    }
    {
      Scoped_Timer fence_timer(timers, "fence");
      L_window.Fence();
      G_window.Fence();
    }
  }
  void Initialize_P_Context::restore_from_residues(
    Vector_Block_Diagonal_Matrix_Residues_Window<double> &bdms_residues,
    std::vector<Block_Diagonal_Matrix> &bdms, Timers &timers,
    El::Matrix<int32_t> &block_timings_ms)
  {
    Scoped_Timer timer(timers, "restore_from_residues");
    for(size_t vec_index = 0; vec_index < bdms.size(); vec_index++)
      {
        auto &bdm = bdms.at(vec_index);
        for(size_t block = 0; block < bdm.blocks.size(); block++)
          {
            auto &matrix = bdm.blocks.at(block);
            const auto &loc = cfg.local_block_locations().at(block);
            const auto node_block_index = loc.block_index_node;
            const auto global_block_index = loc.block_index_global;
            const auto &first_residue_matrix
              = bdms_residues.prime_block_vec_matrices.at(0)
                  .at(node_block_index)
                  .at(vec_index);

            Timer block_timer;
            restore_distmatrix_from_residues(
              std::nullopt, first_residue_matrix, bdms_residues.prime_stride,
              comb, matrix);
            block_timings_ms(global_block_index, 0)
              += block_timer.elapsed_milliseconds();
          }
      }
    Scoped_Timer fence_timer(timers, "fence");
    bdms_residues.Fence();
  }
}