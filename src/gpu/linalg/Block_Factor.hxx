// GPU backend for SDPB: batched Cholesky factors and their explicit inverses
// of equal-size SPD blocks, and P = L^-1 B as products, on device floats.
//
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: the recursion is SDPviBe's chol_inv
// (sdpvibe/linalg.py: L21 = A21 X11^T, S = A22 - L21 L21^T, X21 = -X22 L21 X11)
// and the 32 x 32 leaves are its CUDA-core base case (basecase.py, ported in
// basecase.cuh), each leaf equilibrated by powers of two as in mpf.Chol.
// Portions adapted from SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).
//
// New here: every product uses exact integer accumulation after operand
// alignment (Mpf_Gemm: residues on INT8 tensor cores, then normalization), so the large
// products at the top of the recursion run at tensor-core speed; and
// P = L^-1 B = X B is one batched product, optionally refined once,
// P <- P + X (B - L P).
// Whole-matrix scaling and panel residual correction address the loss of
// backward accuracy observed with the unrefined inverse-product recursion.
#pragma once

#include "../rns/Mpf_Gram.hxx"

#include <map>
#include <memory>
#include <set>
#include <tuple>

namespace sdpb_gpu
{
  // batch blocks of floats, each rows x cols column-major; limb t of entry
  // (i, j) of block b at t*batch*rows*cols + b*rows*cols + j*rows + i
  // A block of each of the batch's matrices can be used in place: then entry
  // (i, j) of block b is at b*bstride + j*ld + i and limb t at t*plane + that
  // (0: packed, ld = rows, bstride = rows*cols, plane = entries()).
  struct Mpf_Batch
  {
    int32_t *m = nullptr;
    int64_t *e = nullptr;
    int L = 0, batch = 0;
    int64_t rows = 0, cols = 0;
    int64_t ld = 0, bstride = 0, plane = 0;
    int64_t entries() const { return int64_t(batch) * rows * cols; }
    int64_t lda() const { return ld ? ld : rows; }
    int64_t bs() const { return bstride ? bstride : rows * cols; }
    int64_t pl() const { return plane ? plane : entries(); }
    Mpf_Matrix view() const
    {
      return Mpf_Matrix{m, e, L, rows, cols, ld, bstride, plane};
    }
    // rows r0..r0+h-1, columns c0..c0+w-1 of every block, in place
    Mpf_Batch block(int64_t r0, int64_t c0, int64_t h, int64_t w) const
    {
      Mpf_Batch B = *this;
      B.m = m + c0 * lda() + r0;
      B.e = e + c0 * lda() + r0;
      B.rows = h;
      B.cols = w;
      B.ld = lda();
      B.bstride = bs();
      B.plane = pl();
      return B;
    }
  };

  struct Factor_Options
  {
    bool equilibrate = true;
    int panel_refine = 1;
    int refine_min_n = 0; // diagnostic: skipping small panels failed real-Q checks
    int correction_digits
      = 0; // 0: full, -1: absolute-error budget, >0: diagnostic
    int symmetric_update_min = 0; // experimental packed Gram update, batch 1
    bool skip_triangular_zeros = true; // exact-zero omission inside residue GEMMs
    bool log_panels = false; // diagnostic: print max exponents of every panel
    // S = A22 - L21 L21^T: lower triangle only, B's residues from A's cache.
    // Nothing above the diagonal of a Schur complement is read later.
    bool lower_update = true;
    // Run-time policy for the panel correction (needs equilibrate). At each
    // panel, the residual that the alignment of A21 and X11 can leave is
    // estimated from their largest exponents,
    //   log2 |A21 - Z L11^T| <~ ceil(log2 n1) + 6 + log2 max|A21| + log2 max|X11| - 8 digits,
    // worst over the batch. While it stays at or below this limit (log2,
    // negative), panels are not corrected. The first panel above it restarts
    // the factorization with every panel corrected (the same result as
    // correction_limit = 0, bit for bit), and later factorizations of this
    // shape correct from the start. 0: every panel is corrected.
    int correction_limit = 0;
    // diagnostic: with log_panels, also form R = A21 - Z L11^T on panels that
    // are not corrected, and log it
    bool measure_residual = false;
    // 0 or the largest panel that is corrected (with correction_limit < 0: as
    // long as the limit holds; after a restart every panel is corrected)
    int refine_max_n = 0;
    // solve(): P = X B by groups of columns of B with the same nonzero rows,
    // found from B's exponents at every call (zero entries of B and of P are
    // not computed; P is bit for bit the same)
    bool structured_solve = true;
    // L21 X11 reads L21^T's residues from L21 L21^T's cache
    bool reuse_residues = true;
    // factor(): crt_guard_bits for its products when > 0 (else the engine's).
    // The truncation is below the worst-case noise of the aligned operands, but
    // not always below the true noise of graded ones: on an inverse used in
    // solves it raised the solver's dual-error floor
    int crt_guard = 0;
    // negative control: the first leaf of this
    // factorization returns its inverse with the low limbs of every mantissa
    // of block 0 cleared, as a leaf would that is accurate only to that many
    // fewer limbs. The MPFR gate must catch it. 0: off.
    int corrupt_leaf_limbs = 0;
  };

  class Block_Factor
  {
  public:
    // batch blocks of size n (n % 32 == 0: pad with identity), L limbs,
    // `digits` per aligned operand in the products
    Block_Factor(int n, int batch, int L, int digits,
                 Gram_Options options = {}, cudaStream_t stream = 0);
    ~Block_Factor();
    // S (n x n, SPD) -> Lf (lower), X = Lf^-1 (lower); upper triangles zero.
    // Precondition, not checked: S is positive definite at the working
    // precision, and so is every Schur complement of the recursion. Unlike
    // El::Cholesky, which throws for a non-HPD matrix, a non-positive pivot
    // is not detected here: the leaves form 1/sqrt of it and the factors are
    // meaningless, with no error.
    // Defaults use global scaling and one full-precision panel correction.
    // Factor_Options{false, 0} reproduces the original unrefined recursion.
    void factor(const Mpf_Batch &S, Mpf_Batch &Lf, Mpf_Batch &X,
                Factor_Options options = {});
    // P = X B (n x ncols), plus `refine` steps P += X (B - Lf P)
    void solve(const Mpf_Batch &Lf, const Mpf_Batch &X, const Mpf_Batch &B,
               Mpf_Batch &P, int refine = 0);
    // C = U V, U^T or V^T optional (allocates C if needed)
    void mm(const Mpf_Batch &U, bool tU, const Mpf_Batch &V, bool tV,
            Mpf_Batch &C, int digits = 0, int triangles = 0);
    // C = U^T V with C's rows of all blocks stacked per column: entry (i, j) of
    // block b at j (batch m) + b m + i, limb planes over batch m n entries
    void mm_flat(const Mpf_Batch &U, const Mpf_Batch &V, int32_t *d_m,
                 int64_t *d_e);
    // the same for U = X^T with X lower triangular, V = B: the entries are bit
    // for bit mm_flat's, by groups of columns of B with the same nonzero rows;
    // what B's and X's zeros make exactly zero is not computed (filled with
    // exact zeros). Costs the scan of B's exponents, one kernel, at every call.
    // guard > 0: crt_guard_bits for this product (then not bit for bit).
    void solve_flat(const Mpf_Batch &Xt, const Mpf_Batch &B, int32_t *d_m,
                    int64_t *d_e, int guard = 0);
    // ms by phase inside the cached products (with options.product_timing)
    void collect_product_timings(std::map<std::string, double> &ms) const;
    Mpf_Batch alloc(int64_t rows, int64_t cols);
    void release(Mpf_Batch &M);
    // with options.timing: milliseconds by stage ("role/...", "copy", "add",
    // "leaf", ...), and with options.profile_products by stage of each product
    // ("kernel/<role>/<stage>"); "count/..." entries count products
    const std::map<std::string, double> &timings() const { return ms_; }
    void reset_timings() { ms_.clear(); }
    size_t cache_evictions() const { return allocation_evictions_; }
    // device memory held by the cached products and their workspace (allocated
    // with cudaMalloc, so not in the stream-ordered pool's high-water mark)
    size_t product_bytes() const
    {
      size_t bytes = opt_.workspace ? opt_.workspace->bytes() : 0;
      for(const auto &g : gemms_)
        bytes += g.second->device_bytes();
      return bytes;
    }
    // correction_limit: panels left uncorrected by the last factor(); restarts
    // since construction; whether shape n x n now corrects from the start
    int panels_skipped() const { return panels_skipped_; }
    int correction_restarts() const { return correction_restarts_; }
    bool corrects_always(int64_t n) const { return corrects_always_.count(n) != 0; }
    // forget those decisions, e.g. at a new call of the solver
    void reset_correction_policy() { corrects_always_.clear(); }
    // free the cached products (they hold buffers sized for their shapes)
    void clear_cache()
    {
      check();
      gemms_.clear();
      grams_.clear();
      bytes_.clear();
      if(workspace_)
        workspace_->release();
    }
    // throws if a product's CRT failed its carry check (synchronizes)
    void check();

  private:
    void rec(const Mpf_Batch &A, const Mpf_Batch &Lf, const Mpf_Batch &X,
             int refine, int min_n, int correction_digits, int gram_min);
    void symmetric_product(const Mpf_Batch &A, Mpf_Batch &C);
    // Lf, X, Z: here and in rec and leaf, parts of the output in place
    void
    refined_panel(const Mpf_Batch &B, const Mpf_Batch &Lf, const Mpf_Batch &X,
                  const Mpf_Batch &Z, int refine, int correction_digits);
    void leaf(const Mpf_Batch &A, const Mpf_Batch &Lf, const Mpf_Batch &X);
    // A'^T B' (mm without its transposes)
    void product(const Mpf_Batch &Ap, const Mpf_Batch &Bp, const Mpf_Batch &C,
                 int digits, int triangles);
    Mpf_Batch transposed(const Mpf_Batch &M);
    void structured_product(const Mpf_Batch &Xt, const Mpf_Batch &B, int32_t *d_m,
                            int64_t *d_e, bool flat);
    // temporaries of factor() and solve(): released buffers are reused for the
    // next allocation of the same size instead of going through the
    // stream-ordered allocator (everything is on one stream, so reuse is
    // ordered like the operations themselves)
    void *pool_get(size_t bytes);
    void pool_put(void *p, size_t bytes);
    void pool_flush();
    Mpf_Gemm &gemm_for(int m, int n, int k, int digits);
    int product_guard_ = 0; // Factor_Options::crt_guard, or solve_flat's guard
    int correction_digits(const Mpf_Batch &R, const Mpf_Batch &X);
    // largest exponents of the nonzero entries of A and B over the batch (as
    // log2 of an upper bound of the magnitudes); LLONG_MIN if all are zero
    void max_tops(const Mpf_Batch &A, const Mpf_Batch &B, long long &ta,
                  long long &tb);
    void tic();
    void toc(const std::string &phase);
    void log_panel(const char *what, const Mpf_Batch &M);
    int n_, batch_, L_, D_;
    // what the planner in gemm_for assumes of a product at this size
    double plan_primes_ = 0, plan_padded_ = 0, plan_crt_digits_ = 0, plan_limbs_ = 0;
    int corrupt_leaf_limbs_ = 0; // Factor_Options::corrupt_leaf_limbs, until used
    bool skip_triangular_zeros_ = false, lower_update_ = false,
         measure_residual_ = false, structured_solve_ = false,
         reuse_residues_ = false;
    int correction_limit_ = 0, panels_skipped_ = 0, correction_restarts_ = 0,
        refine_max_ = 0, refine_min_ = 0;
    std::set<int64_t> corrects_always_;
    std::map<void *, size_t> live_; // temporaries in use during factor()
    size_t allocation_evictions_ = 0;
    Gram_Options opt_;
    cudaStream_t stream_;
    std::map<std::tuple<int, int, int, int, int>, std::unique_ptr<Mpf_Gemm>> gemms_;
    std::map<std::pair<int, int>, std::unique_ptr<Mpf_Gram>> grams_;
    std::map<std::tuple<int, int, int, int, int>, double> bytes_; // estimated buffers
    std::map<std::string, double> ms_;
    std::unique_ptr<Rns_Workspace> workspace_;
    // the leaf's buffers, kept: a leaf is one of the most frequent calls on Q
    uint32_t *leaf_W_ = nullptr, *leaf_Xw_ = nullptr, *leaf_S_ = nullptr,
             *leaf_R_ = nullptr;
    int64_t *leaf_kk_ = nullptr;
    std::multimap<size_t, void *> pool_;
    bool pool_active_ = false;
    std::string role_ = "other";
    bool log_panels_ = false;
    cudaEvent_t ev0_, ev1_;
  };

  // element kernels on Mpf_Batch (all entries, same shapes)
  void mpf_transpose(const Mpf_Batch &A, Mpf_Batch &T, cudaStream_t s);
  // dst = src(r0:r0+h, c0:c0+w) or the reverse (insert = true)
  void mpf_sub_block(const Mpf_Batch &src, int64_t r0, int64_t c0,
                     Mpf_Batch &dst, bool insert, cudaStream_t s);
  // C = A + sgn B (normalized floats)
  void mpf_batch_add(const Mpf_Batch &A, const Mpf_Batch &B, int sgn,
                     const Mpf_Batch &C, cudaStream_t s, bool lower_only = false);
  void mpf_zero(Mpf_Batch &A, cudaStream_t s);
  // M (n x n, symmetric) from its packed lower triangle (packed entries x 1)
  void mpf_unpack_symmetric(const Mpf_Batch &packed, Mpf_Batch &M, cudaStream_t s);
}
