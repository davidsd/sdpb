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
#include <tuple>

namespace sdpb_gpu
{
  // batch blocks of floats, each rows x cols column-major; limb t of entry
  // (i, j) of block b at t*batch*rows*cols + b*rows*cols + j*rows + i
  struct Mpf_Batch
  {
    int32_t *m = nullptr;
    int64_t *e = nullptr;
    int L = 0, batch = 0;
    int64_t rows = 0, cols = 0;
    int64_t entries() const { return int64_t(batch) * rows * cols; }
    Mpf_Matrix view() const { return Mpf_Matrix{m, e, L, rows, cols}; }
  };

  struct Factor_Options
  {
    bool equilibrate = true;
    int panel_refine = 1;
    int refine_min_n = 0; // diagnostic: skipping small panels failed real-Q checks
    int correction_digits
      = 0; // 0: full, -1: absolute-error budget, >0: diagnostic
    int symmetric_update_min = 0; // experimental packed Gram update, batch 1
    bool skip_triangular_zeros = false; // exact-zero omission inside residue GEMMs
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
    Mpf_Batch alloc(int64_t rows, int64_t cols);
    void release(Mpf_Batch &M);
    const std::map<std::string, double> &timings() const { return ms_; }
    size_t cache_evictions() const { return allocation_evictions_; }
    // free the cached products (they hold buffers sized for their shapes)
    void clear_cache()
    {
      check();
      gemms_.clear();
      grams_.clear();
      bytes_.clear();
    }
    // throws if a product's CRT failed its carry check (synchronizes)
    void check();

  private:
    void rec(const Mpf_Batch &A, Mpf_Batch &Lf, Mpf_Batch &X, int refine,
             int min_n, int correction_digits, int gram_min);
    void symmetric_product(const Mpf_Batch &A, Mpf_Batch &C);
    void
    refined_panel(const Mpf_Batch &B, const Mpf_Batch &Lf, const Mpf_Batch &X,
                  Mpf_Batch &Z, int refine, int correction_digits);
    void leaf(const Mpf_Batch &A, Mpf_Batch &Lf, Mpf_Batch &X);
    Mpf_Gemm &gemm_for(int m, int n, int k, int digits);
    int correction_digits(const Mpf_Batch &R, const Mpf_Batch &X);
    void tic();
    void toc(const std::string &phase);
    int n_, batch_, L_, D_;
    bool skip_triangular_zeros_ = false;
    size_t allocation_evictions_ = 0;
    Gram_Options opt_;
    cudaStream_t stream_;
    std::map<std::tuple<int, int, int, int>, std::unique_ptr<Mpf_Gemm>> gemms_;
    std::map<std::pair<int, int>, std::unique_ptr<Mpf_Gram>> grams_;
    std::map<std::tuple<int, int, int, int>, double> bytes_; // estimated buffers
    std::map<std::string, double> ms_;
    cudaEvent_t ev0_, ev1_;
  };

  // element kernels on Mpf_Batch (all entries, same shapes)
  void mpf_transpose(const Mpf_Batch &A, Mpf_Batch &T, cudaStream_t s);
  // dst = src(r0:r0+h, c0:c0+w) or the reverse (insert = true)
  void mpf_sub_block(const Mpf_Batch &src, int64_t r0, int64_t c0,
                     Mpf_Batch &dst, bool insert, cudaStream_t s);
  // C = A + sgn B (normalized floats)
  void mpf_batch_add(const Mpf_Batch &A, const Mpf_Batch &B, int sgn,
                     Mpf_Batch &C, cudaStream_t s);
  void mpf_zero(Mpf_Batch &A, cudaStream_t s);
}
