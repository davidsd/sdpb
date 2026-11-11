// GPU backend for SDPB: SDPB's per-iteration arithmetic on one GPU, below the
// driver SDP_Solver::run.  The driver keeps the
// SIGTERM and checkpoint logic, termination and jump tests, iterations.json
// and the outputs, and sees only scalars.  This header has no CUDA or Elemental
// types, so that SDPB's own compiler can read it.
//
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: blocks of one matrix dimension are batched in
// chunks with inert padding (solver.py Layout, SDP._build), and the passes over
// the chunks follow its Solver.run, which follows SDPB step by step.  Portions
// adapted from SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).
//
// Chunks: the arrays are (batch, rows, cols) device floats in the layout
//   limb t of entry (i, j) of block b at t * batch rows cols + b rows cols + j rows + i
// with, for nb blocks of matrix dimension dim, K = max points, h = max bilinear
// rows, n = dim h, P = K dim (dim + 1) / 2, z = parity nb + q:
//   V (2 nb, h, K), c (nb, P, 1), B (nb, P, N), x (nb, P, 1), X, Y (2 nb, n, n).
// Padding: columns of V beyond the block's K_j and rows beyond its h are zero;
// for each segment s = c (c + 1) / 2 + r of P, the rows s K + k, k >= K_j of c, B,
// x are zero; X and Y are the identity on padded rows.
#pragma once

#include "Host_Float.hxx"
#include "Gpu_Collective.hxx"

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace sdpb_gpu
{
  struct Block_Shape
  {
    int dim, K, he, ho;
  };

  struct Chunk_Info
  {
    std::vector<int> blocks; // global block indices, q = 0 .. nb - 1
    int nb, dim, K, h, P, n;
  };

  enum class Gpu_Array
  {
    V,
    c,
    B,
    x,
    X,
    Y
  };

  struct Gpu_Solver_Options
  {
    Gpu_World world; // default G=1 keeps the original arithmetic path
    int q_tile = 512, q_row_panel = 8192;
    bool q_wide_residual = false;
    int L = 28;
    int digits = 112;
    // relative accuracy of the step lengths (SDPviBe's EIG_TOL)
    double eig_tol = 1e-10;
    // panel residual correction in the factorizations of the blocks of X, Y and
    // S (real Q requires it and always has it)
    bool correct_xy = true, correct_schur = true;
    // kernels of  (0: as before). Adaptive digits of the
    // correction's second product and the correction limit, on Schur blocks; the
    // guard of the truncated CRT in all products; P by column groups
    int schur_correction_digits = 0, schur_correction_limit = 0;
    int crt_guard = 0;
    int pq_guard = 0; // the same guard on P = X B and the factor of Q only
    bool structured_p = false;
    // steps of iterative refinement in the solves with Q (Q stays on the device).
    // One step: p-err 1e-125 -> 4e-158 on real nmax14, CPU SDPB's is 4e-140, at
    // 0.074 s of a 6 s iteration
    int refine_dy = 1;
    // Steps of iterative refinement of the corrector's whole Schur system (S stays
    // on the device for that iteration). Without it the explicit inverses of S leave
    // a floor in the dual residual: 1e-187 on real nmax10, where CPU SDPB reaches
    // 1e-202, so a run with error thresholds of 1e-200 never turns feasible. -1: one
    // step once the dual error is below 2^refine_system_below_log2,
    // so that iterations far from feasibility don't pay for it.
    int refine_system = -1;
    int64_t refine_system_below_log2 = -200;
    // diagnostic, off: refinement of each solve with X in Z and dY (it does not
    // change the floor above)
    int refine_x = 0;
    // at most this many blocks in one chunk (0: one chunk per dimension)
    int max_chunk_blocks = 16;
    // keep each chunk's B as its nonzero entries (80% of B is exact zeros) and
    // expand it on the device for each use; bit for bit the dense B's results
    bool sparse_b = false;
    // release P once Q is formed: the directions use B and the inverse of the
    // Schur factor instead (Zechuan Zheng's lean path); changes the rounding
    bool lean_p = false;
    // per-phase timing inside the exact products (synchronizes)
    bool product_timing = false;
    // deliberate perturbation for the negative control of the comparison with
    // SDPB (relative, 0: none): the corrector's step eigenvalue is multiplied by
    // 1 + perturb_step at iteration perturb_iteration. A double: below about
    // 1e-16 it changes nothing, so it can't test agreement beyond that
    double perturb_step = 0;
    int perturb_iteration = 0;
    // negative control at any depth: at this iteration
    // each chunk's first leaf of the Schur factorization returns its inverse
    // with the low limbs of block 0 cleared (Factor_Options::corrupt_leaf_limbs)
    int corrupt_leaf_limbs = 0;
    int corrupt_leaf_iteration = 0;
  };

  struct Block_Cond
  {
    Host_Float max, min;
  };

  struct Residue_Result
  {
    Host_Float dual_sum;                   // b . y
    std::vector<Host_Float> primal_terms;  // c . x, by chunk
    std::vector<Host_Float> minus_xy_trace; // Tr(-X Y), by chunk
    Host_Float dual_error, primal_error_P, primal_error_p;
  };

  struct Cond_Table
  {
    std::vector<Block_Cond> schur, X, Y; // schur[j], X[2 j + parity], ...
    Block_Cond Q;
  };

  struct Min_Eig
  {
    bool any = false;
    Host_Float value; // the smallest eigenvalue of L^-1 dM L^-T over all blocks
    int refined = 0, unresolved = 0;
  };

  class Gpu_Solver
  {
  public:
    using Refine_Fn = std::function<Host_Float(
      int n, const std::vector<int32_t> &m, const std::vector<int64_t> &e)>;

    Gpu_Solver(const Gpu_Solver_Options &options, int N,
               const std::vector<Block_Shape> &blocks);
    ~Gpu_Solver();

    int num_chunks() const;
    const Chunk_Info &chunk(int i) const;
    // entries of array a of chunk i
    int64_t entries(int i, Gpu_Array a) const;

    // host arrays in the chunk layout above: m[t * entries + idx], e[idx]
    void upload(int i, Gpu_Array a, const int32_t *m, const int64_t *e);
    void download(int i, Gpu_Array a, int32_t *m, int64_t *e);
    void upload_b(const int32_t *m, const int64_t *e);
    void upload_y(const int32_t *m, const int64_t *e);
    void download_y(int32_t *m, int64_t *e);
    // everything uploaded: derived constants, memory plan
    void finish_load();

    // Cholesky factors of X, Y, pairings, residues, objectives, errors, X Y
    Residue_Result residues();
    // Schur blocks and their factors, P = L^-1 B, Q and its factor
    void factor();
    Host_Float R_error(const Host_Float &mu);
    void direction(const Host_Float &beta_mu, bool corrector);
    std::vector<Host_Float> frobenius_terms(); // Tr((X + dX)(Y + dY)), by chunk
    // smallest eigenvalue of L^-1 dM L^-T for the last direction; a block whose
    // float64 eigenvalue with its error bound could move the step length by
    // more than eig_tol goes to refine.
    Min_Eig min_eig(bool dual, double gamma, const Refine_Fn &refine);
    void update(const Host_Float &alpha_P, const Host_Float &alpha_D);
    const Cond_Table &cond() const;
    // CRT checks of all products, profile bookkeeping
    void end_iteration();

    // s by stage: totals and, per iteration, the last iteration (JSON)
    std::string profile_json() const;
    // peak device memory seen by the stages (bytes)
    double peak_device_bytes() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> p_;
  };
}
