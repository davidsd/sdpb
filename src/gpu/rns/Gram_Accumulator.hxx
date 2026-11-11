// GPU backend for SDPB: exact Gram matrix Q = A^T A of a big-integer matrix A
// (rows x N), accumulated over row chunks in residue form, reconstructed once.
//
// This is the GPU counterpart of SDPB's bigint_syrk_blas (src/sdp_solve/
// SDP_Solver/run/bigint_syrk: residues of the normalized P, BLAS per prime, CRT),
// organized as in SDPviBe's GramAcc (sdpvibe/mpf.py: residues of every chunk's
// P^T P accumulated exactly, one CRT per iteration) by Antoine Tilloy
// (https://github.com/AntoineTilloy/SDPviBe), the prototype of this backend.
// Portions adapted from SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).
//
// Input: A as balanced base-256 digits (int8 in [-128, 127]), entry-major,
// entries column-major: digit t of A(k, j) at (j*rows + k)*digits + t.  The
// caller aligns each column to a common exponent first (block floating point),
// as SDPB's Matrix_Normalizer and SDPviBe's GramAcc do.
//
// Pipeline per chunk and group of g primes (all products on INT8 tensor cores,
// cuBLASLt):
//   1. residues: digits (entries x NB) @ T (NB x 2g), T = 256^j mod m_i split
//      hi/lo -> r = 128 z_hi + z_lo mod m_i -> int8 h, l, s = h + l;
//   2. Gram blocks: the lower triangle of Q is cut recursively into off-diagonal
//      rectangles and small diagonal squares; each is one strided-batched INT8
//      GEMM over the 3g residue planes, HH, LL, SS = R_I^T R_J;
//   3. reduce: Q mod m_i += 2^14 HH + 2^7 (SS - HH - LL) + LL  (Karatsuba),
//      kept as uint16 residues of the packed lower triangle.
// finish(): CRT of the packed entries as an INT8 GEMM with the digits of
// M/m_i, then carries -> exact two's complement integers.
#pragma once

#include "Rns_Basis.hxx"

#include <cublasLt.h>
#include <cuda_runtime.h>

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>

namespace sdpb_gpu
{
  class Rns_Workspace;

  struct Gram_Options
  {
    int base = 256;        // diagonal blocks of Q at most this wide (full squares)
    int max_rect = 2048;   // largest Gram GEMM output tile (per side)
    int prime_group = 16;  // primes per residue batch
    int64_t max_entries_per_residue_batch = int64_t(1) << 22;
    int64_t crt_batch = int64_t(1) << 19; // entries per CRT batch
    bool timing = false;   // per-phase CUDA event timing (synchronizes)
    // residues in one fused kernel (WMMA INT8 + modular finish in shared
    // memory) instead of a cuBLASLt GEMM with an int32 intermediate
    bool fused_residues = false;
    // CRT in one fused kernel (WMMA INT8 digit sums in shared memory, carries
    // in-thread) instead of two cuBLASLt GEMMs with int32 intermediates
    bool fused_crt = false;
    // Rns_Gemm: keep the first operand's residues for all primes (3 r bytes
    // per entry) instead of recomputing them for every column panel of B
    bool cache_a = true;
    // Block_Factor: 32 x 32 leaves by the fused Cholesky + inverse kernel
    // (linalg/leaf_fused.cuh) instead of SDPviBe's chol_tiles + trinv_tiles
    bool fused_leaf = true;
    // fused leaf: the next pivot during the update phase (look-ahead)
    bool leaf_lookahead = true;
    // Rns_Gemm: don't wait for the CRT carry check at the end of every product
    // (a host synchronization); the caller checks with Rns_Gemm::check()
    bool defer_crt_check = false;
    // Rns_Gemm: CRT as one INT8 product with base-128 digit tables
    // (Rns_Basis::crt_tables_128) instead of two with base-256 tables and an
    // int32 accumulation pass (the default: exact against GMP, faster)
    bool crt_base128 = true;
    // Mpf_Gemm: round the exact products to floats with the limbs in
    // registers (normalize_blocks_reg) instead of local arrays
    bool fast_rounding = true;
    // Rns_Gemm: the modular reduction of the products also forms the CRT's
    // inputs (digits of y_i and the quotient estimate) in the CRT's row
    // layout, instead of storing uint16 residues that a CRT pre-pass rereads
    bool reduce_to_crt = true;
    // Block_Factor: stage timers inside every cached product as well
    // (synchronizes after each stage; profiles only)
    bool profile_products = false;
    // Block_Factor: stage timers inside its cached products, summed by phase
    // (collect_product_timings)
    bool product_timing = false;
    // Rns_Gemm: when > 0, keep only the top of each exact product (base-128
    // CRT): the digits of C below 2^(8 digits - 2 - guard) are not formed, with
    // an error below 2^(8 digits - 2 - guard) in the integer C. The aligned
    // operands are rounded to integers, so C is uncertain by about k 2^(8 digits - 2)
    // anyway, and a float's unit in C's largest entries is 2^(8 digits - 2):
    // those digits carry nothing. Fewer digits in the CRT product, the carries
    // and the rounding. 0: exact (a negative value cuts into the result: tests).
    int crt_guard_bits = 0;
    // Rns_Gemm: work buffers from this workspace, shared by all products on the
    // stream, instead of each product's own (one product runs at a time)
    Rns_Workspace *workspace = nullptr;
    // Block_Factor: its products share one workspace (the planner then budgets
    // the largest product, not the sum over shapes)
    bool share_workspace = true;
    // Optional rectangular Gram output: columns [0,tile_cols) against
    // [tile_cols,tile_cols+tile_rows). Both lengths must be multiples of 32.
    // Used by distributed Q; zero keeps the original packed triangular path.
    int tile_rows = 0, tile_cols = 0;
  };

  class Gram_Accumulator
  {
  public:
    // N columns; `digits` balanced base-256 digits per entry; rows_total: bound
    // on the total number of rows over all chunks (sets the number of primes).
    Gram_Accumulator(int N, int digits, int64_t rows_total,
                     Gram_Options options = {}, cudaStream_t stream = 0);
    ~Gram_Accumulator();
    Gram_Accumulator(const Gram_Accumulator &) = delete;
    Gram_Accumulator &operator=(const Gram_Accumulator &) = delete;

    // Q += A^T A for a chunk A of `rows` rows (device pointer, layout above).
    void add(const int8_t *d_digits, int64_t rows);

    // Q += A^T A where the caller computes the residue planes itself:
    // fill(i0, gc, R, plane) writes, for primes i0..i0+gc-1 and columns
    // j < N, R[(X*gc + p)*plane + j*rows + k] (X = 0, 1, 2: h, l, s of the
    // centred residue, r = 128 h + l, as in add), zero for padded rows; the
    // padded columns j >= N are zeroed here.  rows % 16 == 0.
    using Residue_Fill
      = std::function<void(int i0, int gc, int8_t *R, int64_t plane)>;
    void add_residues(int64_t rows, const Residue_Fill &fill);
    const int32_t *device_primes() const { return d_primes_; }
    const double *device_minv() const { return d_minv_; }
    cudaStream_t stream() const { return stream_; }

    // Exact lower triangle of Q, packed by columns: entry (i, j), i >= j, has
    // index packed_index(i, j) = j*N - j*(j-1)/2 + (i - j).  d_out: device
    // int32 limb planes, limb t (little endian, two's complement) of entry e at
    // t*packed_size() + e.  Resets the accumulator.
    void finish(int32_t *d_out);
    int64_t packed_size() const
    {
      return opt_.tile_rows ? int64_t(opt_.tile_rows) * opt_.tile_cols
                            : int64_t(N_) * (N_ + 1) / 2;
    }
    // Exact collective seam. Export canonical residues as uint32 (NCCL has
    // no uint16 sum); import reduces the sum modulo each prime. For G ranks
    // require G*(max prime-1) <= UINT32_MAX before calling the collective.
    void export_residues(uint32_t *out) const;
    void import_residues(const uint32_t *sum);
    void reset();
    static int64_t packed_index(int64_t N, int64_t i, int64_t j)
    {
      return j * N - j * (j - 1) / 2 + (i - j);
    }

    int limbs_out() const { return basis_.crt_digits / 4; }
    const Rns_Basis &basis() const { return basis_; }
    int64_t rows_added() const { return rows_added_; }
    // milliseconds by phase (when options.timing)
    const std::map<std::string, double> &timings() const { return ms_; }
    size_t device_bytes() const { return device_bytes_; }
    // INT8 operations issued to the tensor cores so far (2 per multiply-add)
    double int8_ops() const { return ops_; }

  private:
    void gemm_tn(int m, int n, int k, const int8_t *A, int lda, int64_t strideA,
                 const int8_t *B, int ldb, int64_t strideB, int32_t *C, int ldc,
                 int64_t strideC, int batch);
    void gram_rect(int c0, int c1, int i0, int gc, int64_t rows);
    void gram_block(int r0, int r1, int c0, int c1, int i0, int gc, int64_t rows);
    void tic();
    void toc(const std::string &phase);

    int N_, Np_, NB_, g_;
    Rns_Basis basis_;
    Gram_Options opt_;
    cudaStream_t stream_;
    cublasLtHandle_t lt_;
    void *workspace_ = nullptr;
    size_t workspace_bytes_ = size_t(64) << 20;

    // device tables
    int8_t *d_T_ = nullptr;      // [num_padded][2][NB]: hi, lo rows
    int32_t *d_primes_ = nullptr; // [num_padded] (1 in the padding)
    int32_t *d_inv_ = nullptr;   // [num_padded]
    int8_t *d_D_ = nullptr;      // [E][num_padded]
    int32_t *d_Mdig_ = nullptr;  // [E]

    // work buffers
    double *d_minv_ = nullptr;   // [num_padded] 1/m_i
    int32_t *d_resC_ = nullptr;  // residue GEMM output [2g][entries batch]
    int8_t *d_R_ = nullptr;      // residue planes [3][g][Np][rows]
    int64_t rows_cap_ = 0;
    int32_t *d_gram_ = nullptr;  // Gram GEMM output [3g][max_rect^2]
    uint16_t *d_Qres_ = nullptr; // accumulated residues [num primes][packed_size]
    int8_t *d_Yh_ = nullptr, *d_Yl_ = nullptr; // CRT inputs [crt_batch][num_padded]
    int32_t *d_Zh_ = nullptr, *d_Zl_ = nullptr; // CRT digit sums [E][crt_batch]
    int32_t *d_K_ = nullptr;     // CRT quotients [crt_batch]
    int *d_err_ = nullptr;
    int64_t rows_added_ = 0;
    int64_t rows_total_ = 0;
    double ops_ = 0;

    size_t device_bytes_ = 0;
    std::map<std::string, double> ms_;
    cudaEvent_t ev0_, ev1_;
  };
}
