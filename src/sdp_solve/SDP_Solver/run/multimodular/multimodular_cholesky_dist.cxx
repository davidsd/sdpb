#include "multimodular_blas.hxx"
#include "Node_Diag_Exchange.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <optional>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
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
//
// node_mode (--qCholeskyNodeDiag, --qCholeskyNodeFactor; see
// multimodular_blas.hxx and Node_Diag_Exchange.hxx) changes how the
// diagonal blocks reach every rank and how they are factored, not the
// arithmetic: the factor is bit for bit the same.
namespace
{
  using Clock = std::chrono::steady_clock;
  double seconds_since(const Clock::time_point &start)
  {
    return std::chrono::duration<double>(Clock::now() - start).count();
  }

  // --qCholeskyNodeFactor: multimodular_cholesky_lower(A) with its
  // top-level split shared among the node's ranks, bit for bit the same:
  //   A11 := chol(A11)               every rank (BigFloat, n1 <= base_size)
  //   A21^T := A11^{-1} A21^T        node rank r solves the columns c = r mod R
  //                                  (substitution_trsm: columns independent)
  //   A22 -= A21 A21^T (lower)       node rank r updates the rows i = r mod R
  //                                  (multimodular GEMM: each entry is exact
  //                                  before one rounding, with per-row and
  //                                  per-column scales, so any row subset
  //                                  gives the same entries)
  //   A22 := chol(A22)               every rank
  // with the parts shared through node memory (Node_Diag_Exchange). Only
  // where the serial code's triangular solve is substitution_trsm
  // (n/2 <= base_size, e.g. Cholesky_Q's 64-wide panels); otherwise, or on
  // one rank per node, it is multimodular_cholesky_lower itself.
  void node_cholesky_lower(El::Matrix<El::BigFloat> &A, Node_Diag_Exchange &ex,
                           const Multimodular_Blas_Config &config,
                           Multimodular_Blas_Timings &timings)
  {
    const int n = A.Height();
    const int n1 = n / 2, n2 = n - n1;
    const int R = ex.ranks_on_node(), r = ex.rank_on_node();
    if(n <= config.base_size || n1 > config.base_size || R == 1)
      {
        multimodular_cholesky_lower(A, config, timings);
        return;
      }
    const auto top = El::IR(0, n1), bottom = El::IR(n1, n);
    auto A11 = A(top, top);
    auto A21 = A(bottom, top);
    auto A22 = A(bottom, bottom);
    multimodular_cholesky_lower(A11, config, timings);

    // indices q, q + R, q + 2R, ... below m: node rank q's share
    auto share = [R](const int q, const int m) {
      return q < m ? (m - 1 - q) / R + 1 : 0;
    };
    const size_t S = ex.entry_size();
    {
      El::Matrix<El::BigFloat> X(n1, share(r, n2));
      for(int c = r, l = 0; c < n2; c += R, ++l)
        for(int i = 0; i < n1; ++i)
          X(i, l) = A21.CRef(c, i);
      multimodular_trsm_left_lower(A11, X, config, timings);
      const std::byte *buf = ex.node_all_gather(
        [&](const int q) { return size_t(n1) * size_t(share(q, n2)); }, X);
      for(int q = 0; q < R; ++q)
        {
          const auto *p = reinterpret_cast<const El::byte *>(
            buf + ex.node_offset(q) * S);
          for(int c = q; c < n2; c += R)
            for(int i = 0; i < n1; ++i)
              p = A21.Ref(c, i).Deserialize(p);
        }
    }
    {
      const int m = share(r, n2);
      El::Matrix<El::BigFloat> C(m, n2), Ar(m, n1);
      for(int l = 0; l < m; ++l)
        {
          const int i = r + l * R;
          for(int j = 0; j < n2; ++j)
            C(l, j) = A22.CRef(i, j);
          for(int k = 0; k < n1; ++k)
            Ar(l, k) = A21.CRef(i, k);
        }
      // lower triangle: row i >= column j
      std::vector<int> row_start(n2, 0), row_end(n2, -1);
      for(int j = 0; j < n2; ++j)
        {
          const int l0 = j <= r ? 0 : (j - r + R - 1) / R;
          if(l0 < m)
            {
              row_start[j] = l0;
              row_end[j] = m - 1;
            }
        }
      if(m > 0)
        multimodular_gemm_sub(C, Ar, El::TRANSPOSE, A21, row_start, row_end,
                              config, timings);
      const std::byte *buf = ex.node_all_gather(
        [&](const int q) { return size_t(share(q, n2)) * size_t(n2); }, C);
      for(int q = 0; q < R; ++q)
        {
          const auto *p = reinterpret_cast<const El::byte *>(
            buf + ex.node_offset(q) * S);
          for(int j = 0; j < n2; ++j)
            for(int i = q; i < n2; i += R)
              p = A22.Ref(i, j).Deserialize(p);
        }
    }
    multimodular_cholesky_lower(A22, config, timings);
  }
}

void multimodular_cholesky_upper_dist(El::DistMatrix<El::BigFloat> &A,
                                      const int block_size,
                                      const bool mm_diag,
                                      const Multimodular_Blas_Config &config,
                                      Multimodular_Cholesky_Dist_Timings &t,
                                      Timers *timers, const int node_mode)
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
  std::unique_ptr<Node_Diag_Exchange> node_exchange;
  const bool node_diag = (node_mode & mmq_node_diag) != 0,
             node_factor = (node_mode & mmq_node_factor) != 0 && mm_diag;
  if((node_diag || node_factor) && n > 0)
    {
      begin("diag_setup");
      node_exchange = std::make_unique<Node_Diag_Exchange>(
        grid, size_t(Min(bsize, n)) * size_t(Min(bsize, n)));
    }
  std::vector<int> row_start, row_end;
  for(Int k = 0; k < n; k += bsize)
    {
      const Int nb = Min(bsize, n - k);
      const Range<Int> ind1(k, k + nb), ind2(k + nb, n);

      auto A11 = A(ind1, ind1);
      auto A12 = A(ind1, ind2);
      auto A22 = A(ind2, ind2);

      // diag: the copy of the diagonal block to every rank (sub-timer
      // gather, including any wait for the slowest rank's update) and its
      // factorization (factor)
      begin("diag");
      auto start = Clock::now();
      std::optional<Scoped_Timer> sub;
      if(timers != nullptr)
        sub.emplace(*timers, "gather");
      if(node_diag)
        node_exchange->all_gather(A11, A11_STAR_STAR);
      else
        A11_STAR_STAR = A11;
      sub.reset();
      t.diag_gather += seconds_since(start);
      const auto factor_start = Clock::now();
      if(timers != nullptr)
        sub.emplace(*timers, "factor");
      if(mm_diag)
        {
          // L = U^T from the upper triangle; the strictly upper part of
          // L is left as is (multimodular_cholesky_lower reads the lower)
          auto &a = A11_STAR_STAR.Matrix();
          Matrix<BigFloat> L(nb, nb);
          for(Int j = 0; j < nb; ++j)
            for(Int i = 0; i < nb; ++i)
              L(i, j) = i >= j ? a(j, i) : BigFloat(0);
          if(node_factor)
            node_cholesky_lower(L, *node_exchange, config, t.mm);
          else
            multimodular_cholesky_lower(L, config, t.mm);
          for(Int j = 0; j < nb; ++j)
            for(Int i = 0; i <= j; ++i)
              a(i, j) = L(j, i);
        }
      else
        {
          Cholesky(UPPER, A11_STAR_STAR);
        }
      sub.reset();
      t.diag_factor += seconds_since(factor_start);
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
      A12_STAR_MR.AlignWith(A22);
      A12_STAR_MC = A12_STAR_VR;
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
