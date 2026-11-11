#include "Batch_Crt.hxx"

#include "blas_residues.hxx"

#include "residues.hxx"
#include "sdpb_util/assert.hxx"

#include <flint/ulong_extras.h>

#include <cblas.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace
{
  constexpr int chunk_bits = 24;

  // Write sign * (limbs[0..n), little-endian, n significant) as
  // mpf_set_z would store it in an mpf of precision mpf_prec, in
  // El::BigFloat::Serialize's format: int _mp_size, mp_exp_t _mp_exp, then
  // |_mp_size| limbs.
  void serialize_integer(const int sign, const mp_limb_t *limbs, const int n,
                         const int mpf_prec, El::byte *out)
  {
    int size = 0;
    mp_exp_t exp = 0;
    const mp_limb_t *src = limbs;
    if(n > 0)
      {
        const int asize = std::min(n, mpf_prec + 1);
        src = limbs + (n - asize);
        size = sign * asize;
        exp = n;
      }
    std::memcpy(out, &size, sizeof(int));
    out += sizeof(int);
    std::memcpy(out, &exp, sizeof(mp_exp_t));
    out += sizeof(mp_exp_t);
    std::memcpy(out, src, std::abs(size) * sizeof(mp_limb_t));
  }
}

bool Batch_Crt::supported(const Fmpz_Comb &comb)
{
  for(const mp_limb_t p : comb.primes)
    if(!(std::ldexp(double(p) * comb.num_primes, chunk_bits)
         < std::ldexp(1.0, 53)))
      return false;
  return true;
}

Batch_Crt::Batch_Crt(const Fmpz_Comb &comb)
    : K(comb.num_primes), inv_p(K), u(K), primes_d(K), crt_input(K)
{
  for(size_t k = 0; k < K; ++k)
    primes_d[k] = double(comb.primes[k]);
  fmpz_init(z);
  fmpz_t M, M_p;
  fmpz_init_set_ui(M, 1);
  fmpz_init(M_p);
  for(size_t k = 0; k < K; ++k)
    fmpz_mul_ui(M, M, comb.primes[k]);
  const long M_bits = fmpz_bits(M);
  H = (M_bits + chunk_bits - 1) / chunk_bits;
  // sum_t S_t 2^(24 t) < H 2^53 2^(24 (H-1)); q M < K 2^M_bits
  L = (chunk_bits * H + 64 + 8) / GMP_NUMB_BITS + 2;
  M_chunks.assign(K * H, 0);

  mpz_t m;
  mpz_init(m);
  for(size_t k = 0; k < K; ++k)
    {
      const mp_limb_t p = comb.primes[k];
      ASSERT(std::ldexp(double(p) * K, chunk_bits) < std::ldexp(1.0, 53),
             "batch CRT sums would not be exact: p=", p, " num_primes=", K);
      inv_p[k] = 1.0 / double(p);
      fmpz_divexact_ui(M_p, M, p);
      u[k] = n_invmod(fmpz_fdiv_ui(M_p, p), p);
      fmpz_get_mpz(m, M_p);
      for(int t = 0; t < H; ++t)
        {
          mp_limb_t chunk = 0;
          for(int b = 0; b < chunk_bits; ++b)
            chunk |= mp_limb_t(mpz_tstbit(m, chunk_bits * t + b)) << b;
          M_chunks[k + K * t] = double(chunk);
        }
    }
  qM.assign((K + 1) * L, 0);
  for(size_t q = 0; q <= K; ++q)
    {
      fmpz_get_mpz(m, M);
      mpz_mul_ui(m, m, q);
      ASSERT(mpz_size(m) <= size_t(L));
      for(size_t l = 0; l < mpz_size(m); ++l)
        qM[q * L + l] = mpz_getlimbn(m, l);
    }
  acc.resize(L);
  mpz_clear(m);
  fmpz_clear(M_p);
  fmpz_clear(M);
}

Batch_Crt::~Batch_Crt()
{
  fmpz_clear(z);
}

size_t Batch_Crt::restore_serialized(const double *const *residues,
                                     const size_t n, Fmpz_Comb &comb,
                                     const int mpf_prec, El::byte *out,
                                     const size_t slot_size)
{
  if(n == 0)
    return 0;
  ASSERT_EQUAL(comb.num_primes, K);
  Y.resize(n * K);
  S.resize(n * H);
  std::vector<double> quotient(n, 0.0);
  // y_p = (r_p mod p) u_p mod p, and sum_p y_p / p
  for(size_t k = 0; k < K; ++k)
    {
      const nmod_t mod = comb.mods[k];
      const mp_limb_t shift = comb.shifts[k];
      const mp_limb_t uk = u[k];
      const double ip = inv_p[k];
      const double *__restrict r = residues[k];
      double *__restrict y = Y.data() + k * n;
      for(size_t e = 0; e < n; ++e)
        {
          // double_to_uint32_t_residue(r[e], mod), with its shift precomputed
          ASSERT(std::abs(r[e]) <= MAX_BLAS_DP_INT);
          const slong ai = slong(r[e]);
          const mp_limb_t ab = ai < 0 ? mp_limb_t(ai) + shift : mp_limb_t(ai);
          mp_limb_t a, b;
          NMOD_RED(a, ab, mod);
          NMOD_RED(b, a * uk, mod); // a, u < 2^21
          y[e] = double(b);
          quotient[e] += y[e] * ip;
        }
    }
  // S_t, t = 0..H-1, for all n integers
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, H, K, 1.0,
              Y.data(), n, M_chunks.data(), K, 0.0, S.data(), n);

  size_t num_fallback = 0;
  for(size_t e = 0; e < n; ++e)
    {
      El::byte *slot = out + e * slot_size;
      const double qd = quotient[e];
      const double frac = qd - std::floor(qd);
      if(std::abs(frac - 0.5) < fallback_band)
        {
          // ambiguous rounding of q: exact FLINT CRT
          ++num_fallback;
          for(size_t k = 0; k < K; ++k)
            crt_input[k]
              = double_to_uint32_t_residue(residues[k][e], comb.mods[k]);
          fmpz_multi_CRT_ui(z, crt_input.data(), comb.comb, comb.comb_temp,
                            1);
          mpz_t m;
          mpz_init(m);
          fmpz_get_mpz(m, z);
          const int sz = m->_mp_size;
          serialize_integer(sz < 0 ? -1 : 1, m->_mp_d, std::abs(sz),
                            mpf_prec, slot);
          mpz_clear(m);
          continue;
        }
      const size_t q = size_t(std::lround(qd));
      ASSERT(q <= K, DEBUG_STRING(q), DEBUG_STRING(K));
      // c = sum_t S_t 2^(24 t) - q M, two's complement in L limbs
      std::fill(acc.begin(), acc.end(), 0);
      for(int t = 0; t < H; ++t)
        {
          const mp_limb_t v = mp_limb_t(S[e + n * t]);
          const int o = chunk_bits * t;
          int li = o / GMP_NUMB_BITS;
          const int sh = o % GMP_NUMB_BITS;
          const mp_limb_t lo = v << sh;
          const mp_limb_t hi = sh == 0 ? 0 : v >> (GMP_NUMB_BITS - sh);
          mp_limb_t old = acc[li];
          acc[li] += lo;
          mp_limb_t carry = hi + (acc[li] < old);
          while(carry != 0)
            {
              ++li;
              old = acc[li];
              acc[li] += carry;
              carry = acc[li] < old;
            }
        }
      mpn_sub_n(acc.data(), acc.data(), qM.data() + q * L, L);
      int sign = 1;
      if(acc[L - 1] >> (GMP_NUMB_BITS - 1))
        {
          sign = -1;
          mpn_neg(acc.data(), acc.data(), L);
        }
      int size = L;
      while(size > 0 && acc[size - 1] == 0)
        --size;
      serialize_integer(sign, acc.data(), size, mpf_prec, slot);
    }
  return num_fallback;
}

void Batch_Crt::residues(const __mpf_struct *const *x, const size_t n,
                         const long bits, double *out,
                         const size_t row_stride, const size_t prime_stride)
{
  if(n == 0)
    return;
  const int T = (bits + chunk_bits - 1) / chunk_bits;
  if(T != T_res)
    {
      for(size_t k = 0; k < K; ++k)
        ASSERT(std::ldexp(primes_d[k] * T, chunk_bits) < std::ldexp(1.0, 52),
               "batch residues would not be exact: p=", primes_d[k],
               " num_chunks=", T);
      T_res = T;
      W.resize(size_t(T) * K);
      for(size_t k = 0; k < K; ++k)
        {
          const mp_limb_t p = mp_limb_t(primes_d[k]);
          const mp_limb_t base = (mp_limb_t(1) << chunk_bits) % p;
          mp_limb_t power = 1;
          for(int t = 0; t < T; ++t)
            {
              W[t + size_t(T) * k] = double(power);
              power = (power * base) % p;
            }
        }
      chunk_buffer.assign(chunk_bits * T / GMP_NUMB_BITS + 2, 0);
    }
  constexpr size_t batch = 1024;
  chunks.resize(batch * T);
  R.resize(batch * K);
  for(size_t e0 = 0; e0 < n; e0 += batch)
    {
      const size_t nb = std::min(batch, n - e0);
      for(size_t e = 0; e < nb; ++e)
        {
          const __mpf_struct *v = x[e0 + e];
          if(v == nullptr || v->_mp_size == 0)
            for(int t = 0; t < T; ++t)
              chunks[e + t * nb] = 0;
          else
            {
              // |v| < 2^bits, from its top limb
              const long top_bits
                = GMP_NUMB_BITS * v->_mp_exp
                  - __builtin_clzl(v->_mp_d[std::abs(v->_mp_size) - 1]);
              ASSERT(top_bits <= bits, "batch residues: input of ", top_bits,
                     " bits, more than ", bits);
              Blas_Residues::write_chunks(v, 0, T, chunks.data() + e, nb,
                                          chunk_buffer);
            }
        }
      // |sum| < T 2^24 p < 2^52
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nb, K, T, 1.0,
                  chunks.data(), nb, W.data(), T, 0.0, R.data(), nb);
      for(size_t k = 0; k < K; ++k)
        {
          const double p = primes_d[k];
          const double ip = inv_p[k];
          const double half = (p - 1) / 2;
          const double *__restrict r = R.data() + nb * k;
          double *o = out + e0 * row_stride + k * prime_stride;
          for(size_t e = 0; e < nb; ++e)
            {
              // round_mod: [-(p+1)/2, (p+1)/2]; then FLINT's range
              double v = Blas_Residues::round_mod(r[e], p, ip);
              v -= (v > half) ? p : 0.0;
              v += (v < -half) ? p : 0.0;
              o[e * row_stride] = v;
            }
        }
    }
}
