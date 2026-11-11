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

#include <utility>
#include <vector>

namespace sdpb_gpu
{
  // Flags added to Rns_Gemm::multiply_with's `triangles`: the two operands are
  // the same matrix, so B's residues are read from A's cache (needs
  // options.cache_a; m = n, digitsA = digitsB, otherwise ignored); only the
  // lower triangle of C is wanted (the rows of each column panel start at its
  // first column; nothing above is multiplied or written).
  constexpr int Product_Same_Operands = 8, Product_Lower_Only = 16;
  // Mpf_Gemm only: the floats written are -C (the sign is flipped in the
  // rounding, which is symmetric, so this is exactly the negated float)
  constexpr int Product_Negate = 32;

  // One column panel of a structured product (Rns_Gemm::multiply_panels):
  // columns j0..j0+pc-1 of B; rows row0..m-1 of C (those above are known zero);
  // only rows first..last-1 of B are nonzero in these columns (multiples of 16).
  struct Product_Panel
  {
    int j0 = 0, pc = 0, row0 = 0, first = 0, last = 0;
  };

  // Device buffers of a product's work (residues, products, CRT), which matter
  // only during one multiply: products on one stream can share them. Buffers
  // grow to the largest request, which waits for the stream.
  class Rns_Workspace
  {
  public:
    enum Buffer { RA, RA_all, RB, prod, Y_all, Kf, Cres, Y2, Z2, K, exact, lt, count };
    explicit Rns_Workspace(cudaStream_t stream = 0) : stream_(stream) {}
    ~Rns_Workspace();
    Rns_Workspace(const Rns_Workspace &) = delete;
    Rns_Workspace &operator=(const Rns_Workspace &) = delete;
    void *get(Buffer b, size_t bytes);
    size_t bytes() const;
    // gives everything back to the device (waits for the stream)
    void release();

  private:
    cudaStream_t stream_;
    void *p_[count] = {};
    size_t n_[count] = {};
  };

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
    // With Product_Lower_Only, consume() receives rows j0..m-1 of the panel:
    // entry (b, jj, i) at b*(m-j0)*pc + jj*(m-j0) + i - j0, cplane = batch*(m-j0)*pc.
    void multiply_with(const Fill_A &fillA, const Fill_B &fillB, int32_t *d_out,
                       const Consume &consume = nullptr, int triangles = 0);
    // the same over an explicit list of panels: fillB(i0, gc, panel, RB, planeB)
    // writes rows first..last-1 of the panel's pc columns; consume receives rows
    // row0..m-1 (entry (b, jj, i) at b*(m-row0)*pc + jj*(m-row0) + i - row0)
    using Fill_Panel = std::function<void(int i0, int gc, const Product_Panel &,
                                          int8_t *RB, int64_t planeB)>;
    using Consume_Panel = std::function<void(const int32_t *panel, int64_t cplane,
                                             const Product_Panel &)>;
    void multiply_panels(const Fill_A &fillA, const Fill_Panel &fillB,
                         const Consume_Panel &consume,
                         const std::vector<Product_Panel> &panels);
    // for the next multiply only: A's residues are those other's last product
    // cached (the caller guarantees the same A, at the same digits). Ignored
    // if other has no such cache or a different basis or shape.
    void borrow_a_residues(const Rns_Gemm &other);
    // throws if any CRT since the last check failed its carry check
    // (synchronizes; needed only with options.defer_crt_check)
    void check();
    const int32_t *device_primes() const { return d_primes_; }
    const double *device_minv() const { return d_minv_; }
    int panel_width() const { return np_; }
    int m() const { return m_; }
    int n() const { return n_; }
    int k() const { return k_; }
    int batch() const { return batch_; }

    // words of each exact (or, with crt_guard_bits, truncated) result, which
    // is C / 2^exp_shift()
    int limbs_out() const { return Wt_ ? Wt_ : basis_.crt_digits / 4; }
    int exp_shift() const { return 7 * j0_; }
    int64_t plane() const { return int64_t(batch_) * m_ * n_; }
    const Rns_Basis &basis() const { return basis_; }
    double int8_ops() const { return ops_; }
    size_t device_bytes() const { return device_bytes_; }
    const std::map<std::string, double> &timings() const { return ms_; }

  private:
    void residues(const int8_t *d, int64_t ent_per_block, int64_t stride,
                  int digits, const int8_t *T, int i0, int gc, int8_t *R,
                  int64_t plane);
    // points the work buffers into options.workspace
    void bind_workspace();
    // the block products of one prime group are done: reduce them
    void reduce_group(int64_t cplane, int i0, int gc);
    // CRT of the panel in d_prod_'s reductions; the exact result goes to d_out
    // or to consume
    void crt_panel(int64_t cplane, int j0, int pc, int32_t *d_out,
                   const Consume &consume);
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
    const int8_t *borrowed_RA_ = nullptr;
    int32_t *d_prod_ = nullptr;
    uint16_t *d_Cres_ = nullptr;
    int8_t *d_Yh_ = nullptr, *d_Yl_ = nullptr;
    bool reduce_crt_ = false;    // the reduction writes the CRT's rows
    int E2_ = 0;                 // base-128 CRT digits (options.crt_base128)
    int j0_ = 0, Wt_ = 0;        // options.crt_guard_bits: lowest digit kept, words
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
