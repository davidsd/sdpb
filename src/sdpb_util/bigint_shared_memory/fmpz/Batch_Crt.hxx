#pragma once

// Exact CRT of many integers at once, with one dgemm per batch
// Used by restore_and_reduce with
// --batchCrt to restore Q_n from output_residues_window column by
// column, instead of one fmpz_multi_CRT_ui per element.
//
// For an integer c with |c| < M/2 (M = product of the comb's primes) and
// its residues r_p (integer-valued doubles, |r_p| <= 2^53, as the BLAS jobs
// leave them in the window):
//   y_p = (r_p mod p) u_p mod p in [0, p),  u_p = (M/p)^-1 mod p,
//   c = sum_p y_p (M/p) - q M,  q = round(sum_p y_p / p).
// M/p is split into all of its 24-bit chunks m_pt, t = 0..H-1, so that
//   sum_p y_p (M/p) = sum_t S_t 2^(24 t),  S_t = sum_p y_p m_pt,
// exactly: S_t < num_primes 2^21 2^24 < 2^53 is exact in double, and for a
// batch of integers all S_t are one dgemm. (blas_residues.hxx does the same
// with t0 > 0 chunks dropped; here t0 = 0, so the result is exact.)
// q is computed in double: |c|/M < 1/2, and the rounding error of the sum
// is below 1e-12. An element whose sum_p y_p / p lies within 1e-9 of a
// half-integer is restored by FLINT instead (counted in the return value).
//
// The integer is converted to a BigFloat as fmpz_get_mpf does (mpf_set_z:
// the top prec+1 limbs, truncated) and written in El::BigFloat::Serialize's
// format. The result is therefore bit for bit that of
//   restore_bigint_from_residues(); Fmpz_BigInt::to_BigFloat(); Serialize()
// (unit test batch_crt).

#include "Fmpz_Comb.hxx"

#include <El.hpp>

#include <cstddef>
#include <vector>

class Batch_Crt
{
public:
  explicit Batch_Crt(const Fmpz_Comb &comb);

  // The sums S_t are exact only for num_primes p 2^24 < 2^53 for every p.
  // FLINT's primes satisfy it (p < 2^21) when the product needs more than
  // 200 bits, i.e. always for bigint_syrk at SDPB's precisions.
  static bool supported(const Fmpz_Comb &comb);

  // n integers; the residue of integer e modulo prime k is residues[k][e].
  // Writes serialized BigFloats of precision mpf_prec (an mpf's _mp_prec)
  // to out + e * slot_size (slot_size = El::BigFloat::SerializedSize()).
  // Returns the number of elements restored by FLINT (see above).
  size_t restore_serialized(const double *const *residues, size_t n,
                            Fmpz_Comb &comb, int mpf_prec, El::byte *out,
                            size_t slot_size);

  [[nodiscard]] size_t num_primes() const { return K; }

  // --batchResidues (compute_block_residues): residues of the integers
  // c_e = trunc(x_e) (toward zero, as fmpz_set_mpf) modulo the comb's primes,
  // in FLINT's range [-(p-1)/2, (p-1)/2] (that of
  // bigint_to_residues in residues.hxx), so the doubles written are the same.
  // Method: x_e in 24-bit chunks c_et (Blas_Residues::write_chunks), then
  // r_ek = sum_t c_et (2^(24 t) mod p_k) by one dgemm per batch, exact for
  // T 2^24 p < 2^52, T = ceil(bits / 24), and reduced by round_mod.
  // x[e] == nullptr means 0. Requires |x_e| < 2^bits.
  // Output: out[e * row_stride + k * prime_stride].
  void residues(const __mpf_struct *const *x, size_t n, long bits,
                double *out, size_t row_stride, size_t prime_stride);

  // Elements with |frac(sum_p y_p / p) - 1/2| < fallback_band go to FLINT
  // (tests set it to 1 to exercise that path)
  double fallback_band = 1e-9;

private:
  size_t K; // number of primes
  int H;    // number of 24-bit chunks of M (and of M/p)
  int L;    // limbs of the accumulator
  std::vector<double> inv_p;
  std::vector<mp_limb_t> u;
  std::vector<double> M_chunks; // K x H, column-major
  std::vector<mp_limb_t> qM;    // q M, q = 0..K, L limbs each
  std::vector<double> primes_d;
  // residues: W = 2^(24 t) mod p, T_res x K, column-major
  int T_res = 0;
  std::vector<double> W;
  std::vector<double> chunks, R;
  std::vector<mp_limb_t> chunk_buffer;

  // workspace
  std::vector<double> Y, S;
  std::vector<mp_limb_t> acc, crt_input;
  fmpz_t z;

public:
  ~Batch_Crt();
  Batch_Crt(const Batch_Crt &) = delete;
  Batch_Crt &operator=(const Batch_Crt &) = delete;
};
