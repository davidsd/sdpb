#include "multimodular_blas.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <optional>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <vector>

// Cholesky of the global Q (A = U^T U, upper triangle) on its [MC,MR] grid
// This is El::cholesky::UpperVariant3Blocked step for step:
// the same panels, redistributions and panel El::LocalTrsm. The trailing
// update
//
//   A22 -= A12^T A12      (El::LocalTrrk, BigFloat Gemm, ~N^3/(3p) per rank)
//
// is done by multimodular BLAS (multimodular_gemm_sub, see Readme.md) on
// each rank's local part of A22, restricted to the global upper triangle;
// the diagonal blocks too if mm_diag (as sdpb calls it), else by
// El::Cholesky.
// El's BigFloat Gemm costs ~100-200 ns per multiply-add; this is the part
// of Cholesky_Q that grows as N^3/p.
//
// As El, the strictly lower triangle of A is not touched.
namespace
{
  using Clock = std::chrono::steady_clock;
  double seconds_since(const Clock::time_point &start)
  {
    return std::chrono::duration<double>(Clock::now() - start).count();
  }
}

void multimodular_cholesky_upper_dist(El::DistMatrix<El::BigFloat> &A,
                                      const int block_size,
                                      const bool mm_diag,
                                      const Multimodular_Blas_Config &config,
                                      Multimodular_Cholesky_Dist_Timings &t,
                                      Timers *timers)
{
  // Scoped SDPB timers per phase (summed over panels by the profile readers)
  std::optional<Scoped_Timer> phase;
  auto begin = [&](const char *name) {
    phase.reset();
    if(timers != nullptr)
      phase.emplace(*timers, name);
  };
  using namespace El;
  ASSERT_EQUAL(A.Height(), A.Width());
  ASSERT(A.ColDist() == MC && A.RowDist() == MR,
         "multimodular_cholesky_upper_dist: expected an [MC,MR] matrix");
  const Grid &grid = A.Grid();
  DistMatrix<BigFloat, STAR, STAR> A11_STAR_STAR(grid);
  DistMatrix<BigFloat, STAR, VR> A12_STAR_VR(grid);
  DistMatrix<BigFloat, STAR, MC> A12_STAR_MC(grid);
  DistMatrix<BigFloat, STAR, MR> A12_STAR_MR(grid);
  Matrix<BigFloat> A12T;

  const Int n = A.Height();
  const Int bsize = block_size > 0 ? block_size : Blocksize();
  std::vector<int> row_start, row_end;
  for(Int k = 0; k < n; k += bsize)
    {
      const Int nb = Min(bsize, n - k);
      const Range<Int> ind1(k, k + nb), ind2(k + nb, n);

      auto A11 = A(ind1, ind1);
      auto A12 = A(ind1, ind2);
      auto A22 = A(ind2, ind2);

      begin("diag");
      auto start = Clock::now();
      A11_STAR_STAR = A11;
      if(mm_diag)
        {
          // L = U^T from the upper triangle; the strictly upper part of
          // L is left as is (multimodular_cholesky_lower reads the lower)
          auto &a = A11_STAR_STAR.Matrix();
          Matrix<BigFloat> L(nb, nb);
          for(Int j = 0; j < nb; ++j)
            for(Int i = 0; i < nb; ++i)
              L(i, j) = i >= j ? a(j, i) : BigFloat(0);
          multimodular_cholesky_lower(L, config, t.mm);
          for(Int j = 0; j < nb; ++j)
            for(Int i = 0; i <= j; ++i)
              a(i, j) = L(j, i);
        }
      else
        {
          Cholesky(UPPER, A11_STAR_STAR);
        }
      A11 = A11_STAR_STAR;
      t.diag += seconds_since(start);

      begin("panel");
      start = Clock::now();
      A12_STAR_VR.AlignWith(A22);
      A12_STAR_VR = A12;
      LocalTrsm(LEFT, UPPER, ADJOINT, NON_UNIT, BigFloat(1), A11_STAR_STAR,
                A12_STAR_VR);
      t.panel += seconds_since(start);

      begin("redist");
      start = Clock::now();
      A12_STAR_MC.AlignWith(A22);
      A12_STAR_MC = A12_STAR_VR;
      A12_STAR_MR.AlignWith(A22);
      A12_STAR_MR = A12_STAR_VR;
      t.redist += seconds_since(start);

      // Local trailing update, upper triangle in global indices:
      // C(i, j) -= sum_l A12_MC(l, i) A12_MR(l, j) for GlobalRow(i) <=
      // GlobalCol(j). Local rows are increasing in global row, so the rows
      // of local column j are a prefix 0..row_end[j].
      begin("update");
      start = Clock::now();
      {
        auto &C = A22.Matrix();
        const Int m_loc = A22.LocalHeight(), n_loc = A22.LocalWidth();
        if(m_loc > 0 && n_loc > 0)
          {
            const auto &a_mc = A12_STAR_MC.LockedMatrix();
            const auto &a_mr = A12_STAR_MR.LockedMatrix();
            ASSERT_EQUAL(a_mc.Width(), m_loc);
            ASSERT_EQUAL(a_mr.Width(), n_loc);
            A12T.Resize(m_loc, nb);
            for(Int l = 0; l < nb; ++l)
              for(Int i = 0; i < m_loc; ++i)
                A12T(i, l) = a_mc.Get(l, i);
            row_start.assign(n_loc, 0);
            row_end.resize(n_loc);
            Int i = 0;
            for(Int j = 0; j < n_loc; ++j)
              {
                const Int col = A22.GlobalCol(j);
                while(i < m_loc && A22.GlobalRow(i) <= col)
                  ++i;
                row_end[j] = int(i) - 1;
              }
            multimodular_gemm_sub(C, A12T, El::NORMAL, a_mr, row_start,
                                  row_end, config, t.mm);
          }
      }
      t.update += seconds_since(start);

      begin("redist");
      start = Clock::now();
      A12 = A12_STAR_MR;
      t.redist += seconds_since(start);
    }
  phase.reset();
}

size_t
multimodular_cholesky_upper_dist_bytes(const int n, const int grid_height,
                                       const int block_size,
                                       const Multimodular_Blas_Config &config)
{
  if(n <= 0)
    return 0;
  // As in multimodular_block_factor_bytes(): at least 20 bits per prime
  const size_t bits = El::gmp::Precision() + config.guard_bits;
  const size_t num_primes = (2 * bits + 64) / 20 + 1;
  const size_t m_loc = (size_t(n) + grid_height - 1) / std::max(grid_height, 1);
  const size_t nb = block_size > 0 ? block_size : El::Blocksize();
  return sizeof(int32_t) * num_primes * m_loc * (config.column_chunk + nb)
         + sizeof(double) * m_loc * (config.column_chunk + nb);
}
