// GPU backend for SDPB: the Q-formation path on device floats,
//   P (rows x N, per-entry floats)  ->  Q = P^T P (lower triangle, per-entry floats),
// through Gram_Accumulator's exact residue products.
//
// Device floats ("MPF") are SDPviBe's format (sdpvibe/mpf.py, Antoine Tilloy,
// https://github.com/AntoineTilloy/SDPviBe, the prototype of this backend):
// value = m * 2^e with m an L-limb two's complement integer (32-bit limbs,
// little endian, stored as limb planes) normalized to |m| of exactly 32L-2 bits,
// and e an int64 exponent (zero: m = 0, e = MPF_EZ).  Portions adapted from
// SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).
//
// Residues are computed directly from the floats' 32-bit limbs on CUDA cores,
// in one kernel per prime group (align to the column exponent, then
// sum_w a_w (2^(32w) mod m_i) mod m_i), so no digit matrix is materialized.
//
// Accuracy (as SDPB's Matrix_Normalizer and SDPviBe's GramAcc): every entry of
// column j is rounded once to a common exponent, chosen from the column's
// largest entry over all chunks, so that the integers have `digits` balanced
// base-256 digits; the products are exact; each entry of Q is rounded once to
// 32L-2 bits.  The error is normwise per column, |dQ_ij| <~ 2^-(8 digits - 3)
// max_k |P_ki| max_k |P_kj| (rows) .
//
// SDPviBe accumulates chunk by chunk against reference exponents carried from
// the previous iteration, with 28 bits of headroom and a flush when a column
// outgrows them.  Here the chunks of P are resident when Q is formed (the
// backend keeps P), so the exact column maxima are known first and
// no headroom or flush is needed.
#pragma once

#include "Gram_Accumulator.hxx"
#include "Rns_Gemm.hxx"

namespace sdpb_gpu
{
  constexpr int64_t MPF_EZ = -(int64_t(1) << 40);

  // A matrix of device floats, column-major per limb plane (as SDPB's
  // Elemental matrices; SDPviBe's MPF (L, rows, cols) is row-major):
  // m[t*rows*cols + j*rows + i], e[j*rows + i].
  // Mpf_Gemm also takes a part of a larger batch (Mpf_Batch::block): ld is the
  // distance between columns, bstride between blocks, plane between limb planes
  // (0: packed).
  struct Mpf_Matrix
  {
    int32_t *m = nullptr;
    int64_t *e = nullptr;
    int L = 0;
    int64_t rows = 0, cols = 0;
    int64_t ld = 0, bstride = 0, plane = 0;
  };

  class Mpf_Gram
  {
  public:
    // N columns, L limbs per float, `digits` per aligned integer (>= 4L
    // keeps every bit of the column maxima), rows_total over all chunks.
    Mpf_Gram(int N, int L, int digits, int64_t rows_total,
             Gram_Options options = {}, cudaStream_t stream = 0);
    ~Mpf_Gram();
    Mpf_Gram(const Mpf_Gram &) = delete;
    Mpf_Gram &operator=(const Mpf_Gram &) = delete;

    // pass 1: column maxima of every chunk (call for all chunks first)
    void scan(const Mpf_Matrix &P);
    // pass 2: Q += P^T P for each chunk
    void add(const Mpf_Matrix &P);
    // Distributed tile assembly: tops are globally reduced, unsigned-encoded
    // column maxima as returned by mpf_scan_column_tops. map has N entries;
    // -1 denotes a padding column. add_columns also accepts a row-panel view.
    void set_column_tops(const uint64_t *global_tops, const int32_t *map);
    void add_columns(const Mpf_Matrix &P, const int32_t *map);
    // Q's lower triangle, packed by columns (Gram_Accumulator::packed_index),
    // as floats: m[t*packed + e], ex[e].  Resets.
    void finish(int32_t *d_m, int64_t *d_e);

    Gram_Accumulator &gram() { return G_; }
    const std::map<std::string, double> &timings() const { return ms_; }

  private:
    void tic();
    void toc(const std::string &phase);
    int N_, L_, D_, tile_rows_ = 0, tile_cols_ = 0;
    Gram_Accumulator G_;
    cudaStream_t stream_;
    bool timing_;
    int64_t *d_top_ = nullptr;   // column maxima of e + bitlen(m)
    int64_t *d_Eref_ = nullptr;  // column exponents of the aligned integers
    int W_ = 0;                  // table length L + 3
    uint32_t *d_ctab_ = nullptr; // [primes][L+3]: 2^(32 k) mod m_i
    uint32_t *d_ctop_ = nullptr; // [primes][32]: 2^-r mod m_i
    int32_t *d_Qint_ = nullptr;  // exact Q (packed limb planes)
    bool refs_ready_ = false;
    std::map<std::string, double> ms_;
    cudaEvent_t ev0_, ev1_;
  };

  // Caller zeroes N unsigned tops once; scans any number of row chunks and
  // then all-reduces MAX as uint64, never as signed int64.
  void mpf_scan_column_tops(const Mpf_Matrix &P, uint64_t *tops,
                            cudaStream_t stream = 0);

  // Batched exact products of device floats, C_b = A_b^T B_b (A_b: k x m,
  // B_b: k x n), each entry of C rounded once to a float: the same column
  // alignment (per column of A and of B, from exact column maxima), residues
  // from the limbs in registers, Rns_Gemm's block products and CRT.  This is
  // P = X B with an explicit inverse X, and the refinement products.
  // Batched float layout: limb t of entry (r, c) of block b at
  // t*batch*rows*cols + b*rows*cols + c*rows + r (exponents likewise without t).
  class Mpf_Gemm
  {
  public:
    Mpf_Gemm(int m, int n, int k, int batch, int L, int digits,
             Gram_Options options = {}, int panel = 512,
             cudaStream_t stream = 0);
    ~Mpf_Gemm();
    Mpf_Gemm(const Mpf_Gemm &) = delete;
    Mpf_Gemm &operator=(const Mpf_Gemm &) = delete;
    // A: rows k, cols m; B: rows k, cols n (batch blocks each); C: m x n floats
    // triangles: see Rns_Gemm::multiply_with. With Product_Lower_Only the
    // entries of C above the diagonal are not written (panel by panel: rows
    // above a panel's first column; those of the panel's own triangle are).
    // flat: C's rows of all blocks stacked per column (the layout Mpf_Gram
    // takes), entry (i, j) of block b at j (batch m) + b m + i
    void multiply(const Mpf_Matrix &A, const Mpf_Matrix &B, int32_t *d_cm,
                  int64_t *d_ce, int triangles = 0, bool flat = false);
    // C = A^T B with A's columns triangular (A^T lower) and B known to be
    // zero except on groups of columns sharing their nonzero rows: for each
    // group, cols are its columns of B, rows first..last-1 (multiples of 16)
    // hold all its nonzero entries. Rows of C above `first` are known zero and
    // are not written. Bit for bit the entries of multiply().
    struct Column_Group
    {
      std::vector<int> cols;
      int first = 0, last = 0;
    };
    void multiply_groups(const Mpf_Matrix &A, const Mpf_Matrix &B,
                         const std::vector<Column_Group> &groups, int32_t *d_cm,
                         int64_t *d_ce, bool flat = false);
    // for the next multiply: C is part of a larger batch, with these distances
    // between limb planes, blocks and columns (default: packed)
    void set_output_layout(int64_t plane, int64_t bstride, int64_t ld)
    {
      out_plane_ = plane;
      out_bstride_ = bstride;
      out_ld_ = ld;
    }
    Rns_Gemm &gemm() { return G_; }
    const Rns_Gemm &gemm() const { return G_; }
    const std::map<std::string, double> &timings() const { return ms_; }
    size_t device_bytes() const { return G_.device_bytes(); }

  private:
    // column maxima of A and B and the alignment exponents
    void alignments(const Mpf_Matrix &A, const Mpf_Matrix &B, bool same);
    void tic();
    void toc(const std::string &phase);
    int m_, n_, k_, kp_, mp_, batch_, L_, D_;
    bool fast_rounding_ = true;
    Rns_Gemm G_;
    cudaStream_t stream_;
    bool timing_;
    uint32_t *d_ctab_ = nullptr;
    int64_t *d_topA_ = nullptr, *d_topB_ = nullptr, *d_EA_ = nullptr,
            *d_EB_ = nullptr;
    int32_t *d_exact_ = nullptr;
    int32_t *d_colmap_ = nullptr; // position in the groups -> column of B
    int64_t out_plane_ = 0, out_bstride_ = 0, out_ld_ = 0;
    std::map<std::string, double> ms_;
    cudaEvent_t ev0_, ev1_;
  };
}
