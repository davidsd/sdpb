#pragma once

// Blocked Cholesky and triangular solve for El::Matrix<El::BigFloat> whose
// O(n^3) work is done as matrix multiplications on fixed-point integers:
// residues modulo ~100 primes, cblas_dgemm for each prime, and CRT,
// as in bigint_syrk (see bigint_syrk/Readme.md); residues and CRT by BLAS
// (sdpb_util/bigint_shared_memory/fmpz/blas_residues.hxx). See Readme.md.
// Small diagonal blocks are handled in BigFloat.
//
// Precision: in each product C -= A op(B), row i of A and column j of op(B)
// are scaled by powers of 2 so that their largest entries are ~2^bits, with
// bits = El::gmp::Precision() + guard_bits, and rounded to integers.
// The integer product is exact, so the error of C_ij is bounded by
// ~k * 2^-bits * max_l |A_il| * max_l |op(B)_lj|.
// For Cholesky this is the usual normwise backward error bound
// |dA_ij| <~ eps * sqrt(A_ii A_jj).

#include <El.hpp>

#include <cstddef>
#include <vector>

// Test hook: SDPB_MULTIMODULAR_BASE_SIZE and SDPB_MULTIMODULAR_DISTRIBUTED_MIN
// (environment) override base_size and distributed_min_size, so that the
// multimodular kernels run on small test SDPs, whose Schur blocks are below
// the default base_size.
struct Multimodular_Blas_Config;
void apply_multimodular_test_overrides(Multimodular_Blas_Config &config);

// Throws unless the GEMM and Trsm residues and CRT are exact at this
// precision (bits), for any inner dimension up to 2^24: up to 3200 bits
// with the default config. sdpb checks it at startup.
void check_multimodular_precision(const Multimodular_Blas_Config &config,
                                  long precision);

struct Multimodular_Blas_Config
{
  // Diagonal blocks of at most this size are processed in BigFloat
  // (Cholesky / substitution); also the panel height of the
  // left-looking Trsm
  int base_size = 32;
  // Multimodular_Block_Factor: Cholesky levels of at least this many rows
  // are split over the block's grid
  int distributed_min_size = 256;
  // Extra bits of the fixed-point integers on top of El::gmp::Precision()
  int guard_bits = 64;
  // Columns of op(B) processed at once (bounds the memory for residues)
  int column_chunk = 256;
  // If > 0: the residues of A in C -= A op(B) may use at most this memory;
  // larger A are processed in row blocks (op(B) is then converted once per
  // block)
  size_t gemm_A_bytes = 0;
  // Reconstruct products with FLINT's fmpz_multi_CRT_ui (exact, slower)
  // instead of the truncated BLAS CRT. For testing.
  bool flint_crt = false;
  // Trsm: left-looking (each conversion once) instead of recursive, when
  // B has at most left_looking_max_ratio * L.Height() columns. (Measured:
  // left-looking is faster for few columns, e.g. a 975-row block with 262
  // columns per rank; recursive is faster for many, e.g. 195 x 1834.)
  bool left_looking = true;
  double left_looking_max_ratio = 2.0;
  // Left-looking Trsm: headroom of the per-column fixed-point scale of X
  int growth_bits = 32;
  // Left-looking Trsm: memory for the stored residues of X; sets the
  // number of columns processed at once (L is converted once per chunk)
  size_t left_looking_bytes = size_t(384) << 20;
};

struct Multimodular_Blas_Timings
{
  double base_seconds = 0;     // El::Cholesky / El::Trsm on diagonal blocks
  double to_int_seconds = 0;   // scaling, BigFloat -> fmpz, residues
  double blas_seconds = 0;     // cblas_dgemm over all primes
  double from_int_seconds = 0; // CRT, fmpz -> BigFloat, update of C
  size_t num_gemms = 0;
  size_t max_num_primes = 0;
  size_t num_rescales = 0; // left-looking Trsm: column scale raised
  void add(const Multimodular_Blas_Timings &other);
};

// C -= A * op(B), op(B) = B or B^T.
// If lower_only, only C(i,j) with i >= j are updated
// (the strictly upper triangle of C is left untouched).
void multimodular_gemm_sub(El::Matrix<El::BigFloat> &C,
                           const El::Matrix<El::BigFloat> &A,
                           El::Orientation orientation_B,
                           const El::Matrix<El::BigFloat> &B, bool lower_only,
                           const Multimodular_Blas_Config &config,
                           Multimodular_Blas_Timings &timings);
// Same, but only C(i,j) with i >= row_start[j] are updated.
void multimodular_gemm_sub(El::Matrix<El::BigFloat> &C,
                           const El::Matrix<El::BigFloat> &A,
                           El::Orientation orientation_B,
                           const El::Matrix<El::BigFloat> &B,
                           const std::vector<int> &row_start,
                           const Multimodular_Blas_Config &config,
                           Multimodular_Blas_Timings &timings);

// B := L^{-1} B, L lower triangular with nonzero diagonal
// (only the lower triangle of L is read).
void multimodular_trsm_left_lower(const El::Matrix<El::BigFloat> &L,
                                  El::Matrix<El::BigFloat> &B,
                                  const Multimodular_Blas_Config &config,
                                  Multimodular_Blas_Timings &timings);

// Left-looking variant (see multimodular_gemm.cxx); used by
// multimodular_trsm_left_lower when config.left_looking.
void multimodular_trsm_left_looking(const El::Matrix<El::BigFloat> &L,
                                    El::Matrix<El::BigFloat> &B,
                                    const Multimodular_Blas_Config &config,
                                    Multimodular_Blas_Timings &timings);

// B := L^{-1} B by forward substitution in BigFloat (base case)
void substitution_trsm(const El::Matrix<El::BigFloat> &L,
                       El::Matrix<El::BigFloat> &B);

// Lower Cholesky factor in BigFloat for small A, lower triangle only
void base_cholesky_lower(El::Matrix<El::BigFloat> &A);

// A := L, the lower Cholesky factor of A, as El::Cholesky(LOWER, A):
// only the lower triangle is read and written.
void multimodular_cholesky_lower(El::Matrix<El::BigFloat> &A,
                                 const Multimodular_Blas_Config &config,
                                 Multimodular_Blas_Timings &timings);
