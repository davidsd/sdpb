// GPU backend for SDPB: exact batched products C_b = A_b^T B_b of big-integer
// matrices of equal shape, in residue arithmetic on INT8 tensor cores.
//
// This is the engine for SDPB's per-block products (L_j^-1 B_j applied as a
// product with an explicit inverse, bilinear pairings, refinement steps),
// batched over blocks of one size as in SDPviBe's batched.py (bgemm_nt) by
// Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe), the prototype of
// this backend.  Portions adapted from SDPviBe, Copyright (c) 2026 Antoine
// Tilloy, MIT License (docs/GPU_backend_credits.md).
//
// Layout (as Gram_Accumulator): a (rows x cols) block is stored column-major
// by entries, each entry `digits` balanced base-256 digits: digit t of entry
// (r, c) of block b at ((b*block_stride + c*rows + r)*digits + t).
// A_b is k x m, B_b is k x n; C_b (m x n) is written exactly as two's
// complement 32-bit limb planes: limb w of C_b(i, j) at
// w*plane + b*m*n + j*m + i, plane = batch*m*n.
//
// New here: three INT8 products per prime (Karatsuba), column panels of B so
// that the residues of C stay bounded, one CRT per panel.
#pragma once

#include "Gram_Accumulator.hxx"

namespace sdpb_gpu
{
  class Rns_Gemm
  {
  public:
    Rns_Gemm(int m, int n, int k, int batch, int digitsA, int digitsB,
             Gram_Options options = {}, int panel = 512,
             cudaStream_t stream = 0);
    ~Rns_Gemm();
    Rns_Gemm(const Rns_Gemm &) = delete;
    Rns_Gemm &operator=(const Rns_Gemm &) = delete;

    // block strides in entries (default: k*m for A, k*n for B)
    void multiply(const int8_t *dA, const int8_t *dB, int32_t *d_out,
                  int64_t strideA = -1, int64_t strideB = -1);

    // the same with caller-computed residue planes: fillA(i0, gc, RA, planeA)
    // writes, for primes i0..i0+gc-1, the h, l, s planes
    // RA[(X*gc + p)*planeA + b*k*m + i*k + t] of every block's A;
    // fillB(i0, gc, j0, pc, RB, planeB) those of columns j0..j0+pc-1 of B,
    // RB[(X*gc + p)*planeB + b*k*pc + jj*k + t].  With options.cache_a, A's
    // residues are computed once for all panels.
    using Fill_A = std::function<void(int i0, int gc, int8_t *RA, int64_t planeA)>;
    using Fill_B = std::function<void(int i0, int gc, int j0, int pc, int8_t *RB,
                                      int64_t planeB)>;
    // d_out == nullptr: each panel's exact result (W limb planes over
    // batch*m*pc entries, entry (b, jj, i) at b*m*pc + jj*m + i) is handed to
    // consume(panel, cplane, j0, pc) instead, and nothing of size m*n is kept
    using Consume = std::function<void(const int32_t *panel, int64_t cplane,
                                       int j0, int pc)>;
    // triangles: 0 general, 1 logical left factor A^T lower triangular,
    // 2 B lower triangular, 4 B upper triangular. The caller guarantees the
    // omitted entries are exactly zero. Alignment and CRT are unchanged.
    void multiply_with(const Fill_A &fillA, const Fill_B &fillB, int32_t *d_out,
                       const Consume &consume = nullptr, int triangles = 0);
    // throws if any CRT since the last check failed its carry check
    // (synchronizes; needed only with options.defer_crt_check)
    void check();
    const int32_t *device_primes() const { return d_primes_; }
    const double *device_minv() const { return d_minv_; }
    int m() const { return m_; }
    int n() const { return n_; }
    int k() const { return k_; }
    int batch() const { return batch_; }

    int limbs_out() const { return basis_.crt_digits / 4; }
    int64_t plane() const { return int64_t(batch_) * m_ * n_; }
    const Rns_Basis &basis() const { return basis_; }
    double int8_ops() const { return ops_; }
    size_t device_bytes() const { return device_bytes_; }
    const std::map<std::string, double> &timings() const { return ms_; }

  private:
    void residues(const int8_t *d, int64_t ent_per_block, int64_t stride,
                  int digits, const int8_t *T, int i0, int gc, int8_t *R,
                  int64_t plane);
    void tic();
    void toc(const std::string &phase);

    int m_, n_, k_, batch_, DA_, DB_, g_, np_;
    Rns_Basis basis_;
    Gram_Options opt_;
    cudaStream_t stream_;
    cublasLtHandle_t lt_;
    void *workspace_ = nullptr;
    size_t workspace_bytes_ = size_t(64) << 20;
    int8_t *d_TA_ = nullptr, *d_TB_ = nullptr;
    int32_t *d_primes_ = nullptr, *d_inv_ = nullptr, *d_Mdig_ = nullptr;
    double *d_minv_ = nullptr;
    int8_t *d_D_ = nullptr;
    int32_t *d_resC_ = nullptr;
    int8_t *d_RA_ = nullptr, *d_RB_ = nullptr;
    int8_t *d_RAall_ = nullptr; // cached residues of A, all primes
    int32_t *d_prod_ = nullptr;
    uint16_t *d_Cres_ = nullptr;
    int8_t *d_Yh_ = nullptr, *d_Yl_ = nullptr;
    int E2_ = 0;                 // base-128 CRT digits (options.crt_base128)
    int8_t *d_D2_ = nullptr, *d_Y2_ = nullptr;
    int32_t *d_Mdig2_ = nullptr, *d_Z2_ = nullptr;
    int8_t *d_Yall_ = nullptr;    // CRT rows of a panel (options.reduce_to_crt)
    float *d_Kf_ = nullptr, *d_minvf_ = nullptr;
    int32_t *d_Zh_ = nullptr, *d_Zl_ = nullptr, *d_K_ = nullptr;
    int *d_err_ = nullptr;
    double ops_ = 0;
    size_t device_bytes_ = 0;
    std::map<std::string, double> ms_;
    cudaEvent_t ev0_, ev1_;
  };
}
