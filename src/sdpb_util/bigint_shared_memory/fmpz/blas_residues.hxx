#pragma once

// Residues and CRT computed with BLAS, for fixed-point integers taken
// directly from BigFloat (mpf) values, without intermediate fmpz.
// They complement the per-element routines of residues.hxx (FLINT's
// fmpz_multi_mod_ui, fmpz_multi_CRT_ui). The residues are the same as FLINT's. The CRT is
// exact only for small products (t0 = 0 below); otherwise its result has
// the bounded error stated below:
//
// Residues. An integer x = trunc(|v| 2^shift) sign(v), |x| < 2^bits, is
//   split into 24-bit chunks c_t, and
//     x mod p = sum_t c_t (2^(24t) mod p)  mod p,
//   which is exact in double for T 2^24 p < 2^53 (T = number of chunks).
//   For a batch of integers this is one dgemm (chunks x powers of 2^24).
//
// CRT. With primes p_1..p_n of product M, |c| < M/8, and
//   y_p = c u_p mod p in [0, p), u_p = (M/p)^-1 mod p:
//     c = sum_p y_p (M/p) - q M,  q = round(sum_p y_p / p),
//   where q is computed in double (|c|/M < 1/8 leaves room for its rounding
//   error). M/p is split into 24-bit chunks m_pt; S_t = sum_p y_p m_pt is
//   exact in double, and one dgemm gives all S_t for a batch of integers.
//   Only chunks t >= t0 are formed (t0 = 0, i.e. exact, for
//   min(bits_A, bits_B) < 141): the result c' (in units of 2^(24 t0))
//   satisfies
//   -2^(24 t0) < c - c' 2^(24 t0) < num_primes 2^21 2^(24 t0),
//   with 24 t0 <= min(bits_A, bits_B) - 117. For a product of fixed-point
//   integers |a| < 2^bits_A, |b| < 2^bits_B this is at least 2^89 below the
//   error of truncating the factors to integers (up to 2^bits_B per term).
//   Conversion to an mpf of precision P (blas_crt_to_mpf) truncates
//   c' 2^(24 t0) once, so the result v satisfies
//     |v - c| < 2 num_primes 2^21 2^(24 t0) + 2^-P |c|.
//   This bound is the contract (unit test blas_residues). The result is
//   not bitwise equal to the conversion of the exact c (flint_crt_to_mpf)
//   in general. t0 is a multiple of 8, so that 24 t0 is a multiple of the
//   limb size and both conversions truncate at the same bit; they usually
//   agree, but differ where c lies within the CRT error of a truncation
//   boundary:
//   - the error borrows across the kept limbs when c has a long run of zero
//     bits there; e.g. for c = 2^K, v is one unit of the last kept limb
//     below FLINT's (at P = 768: 2^-894 relative for K = 1662, and in
//     multimodular_gemm_sub, 1 - 1*1 gives 2^-830 instead of 0);
//   - |c| is so small (< ~2^(24 t0 + 28 + P)) that the error reaches the
//     kept limbs.
//   Both are within the bound above, i.e. below the working precision.
//   The multiplication by u_p is folded into the residues of one factor
//   (multiply_by_u in compute_blas_residues).
//
// Used by the multimodular Cholesky/Trsm (sdp_solve/SDP_Solver/run/multimodular).

#include "Fmpz_Comb.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <cblas.h>

#include <cstdint>
#include <vector>

namespace Blas_Residues
{
  constexpr int chunk_bits = 24;

  // v - round(v/p) p, in [-(p+1)/2, (p+1)/2], for an integer-valued double
  // |v| < 2^52. Branch-free (rounding by the 1.5*2^52 trick) so that the
  // loops vectorize.
  inline double round_mod(const double v, const double p, const double inv_p)
  {
    constexpr double magic = 6755399441055744.0; // 1.5 * 2^52
    const double q = (v * inv_p + magic) - magic;
    return v - q * p;
  }
  // v mod p in [0, p)
  inline double mod_p(const double v, const double p, const double inv_p)
  {
    const double r = round_mod(v, p, inv_p);
    return r + ((r < 0) ? p : 0.0);
  }

  // Chunks [chunk_bits t, chunk_bits (t+1)) of the integer
  // y = trunc(|x| * 2^shift) < 2^(chunk_bits num_chunks), times the sign of
  // x, for t = 0..num_chunks-1 (written with stride). buffer: scratch space
  // of at least chunk_bits * num_chunks / 64 + 2 limbs.
  void
  write_chunks(const __mpf_struct *x, long shift, int num_chunks,
               double *chunks, size_t stride, std::vector<mp_limb_t> &buffer);

  // Primes and BLAS tables for products of integers |a| < 2^bits_A,
  // |b| < 2^bits_B with inner dimension <= k.
  // Requires bits_A + bits_B + FLINT_BIT_COUNT(k) + 3 > 200, so that
  // Fmpz_Comb's primes are below 2^21 (exactness of the chunked sums).
  // Product bound: |c| < k 2^(bits_A + bits_B); 3 more bits make
  // |c|/M < 1/8 for the rounding of q.
  struct Context
  {
    Fmpz_Comb comb; // the primes; FLINT CRT (flint_crt_to_mpf)
    size_t num_primes;
    std::vector<double> primes, inverse_primes;
    std::vector<double> u; // (M/p)^-1 mod p
    // residues: |x| < 2^max(bits_A, bits_B) in num_chunks chunks;
    // W = 2^(24t) mod p, num_chunks x num_primes, column-major
    int num_chunks;
    std::vector<double> W;
    // CRT: chunks t0..t0+num_high_chunks-1 of M/p,
    // num_primes x num_high_chunks, column-major
    int t0, num_high_chunks;
    std::vector<double> M_chunks;
    // floor(q M / 2^(24 t0)), q = 0..num_primes, num_limbs limbs each
    int num_limbs;
    std::vector<mp_limb_t> qM;

    Context(long bits_A, long bits_B, slong k);
  };

  // A Context for each (bits_A, bits_B, k), created once per process
  const Context &get_context(long bits_A, long bits_B, slong k);

  // Whether the chunked residues and CRT of Context(bits_A, bits_B, k) are
  // exact in double: num_chunks * p * 2^24 < 2^53 and
  // num_primes * p * 2^24 < 2^53 for every prime p. The Context constructor
  // fails otherwise. With the default Multimodular_Blas_Config this limits
  // --multimodularTrsm to --precision=3200 (check_multimodular_precision).
  bool is_exact(long bits_A, long bits_B, slong k);

  // Residues of n fixed-point integers trunc(x_e * 2^shift_e), stored as
  // int32: out[e + prime_stride * i] = residue modulo prime i in
  // [-(p_i+1)/2, (p_i+1)/2], multiplied by u_i if multiply_by_u.
  // element(e, shift) returns x_e (nullptr for zero) and sets shift.
  template <class Element>
  void compute_blas_residues(const Context &context, size_t n,
                             const Element &element, bool multiply_by_u,
                             int32_t *out, size_t prime_stride = 0)
  {
    if(prime_stride == 0)
      prime_stride = n;
    const size_t num_primes = context.num_primes;
    const int T = context.num_chunks;
    constexpr size_t batch = 2048;
    std::vector<double> chunks(batch * T), residues(batch * num_primes);
    std::vector<mp_limb_t> buffer(chunk_bits * T / GMP_NUMB_BITS + 2);
    for(size_t e0 = 0; e0 < n; e0 += batch)
      {
        const size_t nb = std::min(batch, n - e0);
        for(size_t e = 0; e < nb; ++e)
          {
            long shift;
            const __mpf_struct *x = element(e0 + e, shift);
            if(x == nullptr || x->_mp_size == 0)
              {
                for(int t = 0; t < T; ++t)
                  chunks[e + t * nb] = 0;
              }
            else
              write_chunks(x, shift, T, chunks.data() + e, nb, buffer);
          }
        // |sum| < T 2^24 p < 2^53
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nb, num_primes,
                    T, 1.0, chunks.data(), nb, context.W.data(), T, 0.0,
                    residues.data(), nb);
        for(size_t i = 0; i < num_primes; ++i)
          {
            const double p = context.primes[i];
            const double inv_p = context.inverse_primes[i];
            const double *__restrict r = residues.data() + nb * i;
            int32_t *__restrict o = out + e0 + prime_stride * i;
            if(multiply_by_u)
              {
                const double u = context.u[i];
                for(size_t e = 0; e < nb; ++e)
                  o[e] = int32_t(
                    round_mod(round_mod(r[e], p, inv_p) * u, p, inv_p));
              }
            else
              {
                for(size_t e = 0; e < nb; ++e)
                  o[e] = int32_t(round_mod(r[e], p, inv_p));
              }
          }
      }
  }

  // Scratch space for the CRT
  struct Crt_Workspace
  {
    std::vector<double> Y, S;
    std::vector<mp_limb_t> crt_input, acc;
    fmpz_t z;
    explicit Crt_Workspace(const Context &context);
    ~Crt_Workspace();
    Crt_Workspace(const Crt_Workspace &) = delete;
    Crt_Workspace &operator=(const Crt_Workspace &) = delete;
  };

  // BLAS CRT, step 1: for n integers with residues y_p (in [0, p), already
  // multiplied by u_p) at residues[e + residue_stride * p], compute
  // work.S[e + n t] = S_t, t = 0..num_high_chunks-1 (one dgemm).
  void
  blas_crt_high_chunks(const Context &context, const int32_t *residues,
                       size_t residue_stride, size_t n, Crt_Workspace &work);

  // BLAS CRT, step 2, for one integer: value := c' 2^(24 t0) 2^shift,
  // truncated to the precision P of value, i.e. c 2^shift within
  // (2 num_primes 2^21 2^(24 t0) + 2^-P |c|) 2^shift (see above),
  // with S_t = S[t * stride] from blas_crt_high_chunks
  // and quotient = sum_p y_p / p. Returns false (value untouched) if the
  // result is zero.
  bool blas_crt_to_mpf(const Context &context, const double *S, size_t stride,
                       double quotient, long shift, Crt_Workspace &work,
                       mpf_t value);

  // Exact CRT with FLINT's fmpz_multi_CRT_ui (residues y_p in [0, p), not
  // multiplied by u_p, at residues[residue_stride * p]):
  // value := c 2^shift. Returns false (value untouched) if c = 0.
  bool flint_crt_to_mpf(const Context &context, const int32_t *residues,
                        size_t residue_stride, long shift, Crt_Workspace &work,
                        mpf_t value);
}
