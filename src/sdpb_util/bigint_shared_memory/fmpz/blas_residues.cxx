#include "blas_residues.hxx"

#include <flint/ulong_extras.h>

#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <tuple>

namespace Blas_Residues
{
  void write_chunks(const __mpf_struct *x, const long shift,
                    const int num_chunks, double *chunks, const size_t stride,
                    std::vector<mp_limb_t> &buffer)
  {
    const long size = std::abs(x->_mp_size);
    const mp_limb_t *d = x->_mp_d;
    const double sign = x->_mp_size < 0 ? -1.0 : 1.0;
    const long num_limbs = long(chunk_bits) * num_chunks / GMP_NUMB_BITS + 2;
    std::fill(buffer.begin(), buffer.begin() + num_limbs, 0);
    // bit position in y of bit 0 of d[0]
    const long offset = long(GMP_NUMB_BITS) * (x->_mp_exp - size) + shift;
    if(offset >= 0)
      {
        const long ls = offset / GMP_NUMB_BITS;
        const unsigned bs = offset % GMP_NUMB_BITS;
        ASSERT(ls + size + 1 <= num_limbs, "fixed-point overflow");
        if(bs == 0)
          std::copy(d, d + size, buffer.data() + ls);
        else
          buffer[ls + size] = mpn_lshift(buffer.data() + ls, d, size, bs);
      }
    else
      {
        const long ld = (-offset) / GMP_NUMB_BITS;
        const unsigned bs = (-offset) % GMP_NUMB_BITS;
        if(ld < size)
          {
            ASSERT(size - ld <= num_limbs, "fixed-point overflow");
            if(bs == 0)
              std::copy(d + ld, d + size, buffer.data());
            else
              mpn_rshift(buffer.data(), d + ld, size - ld, bs);
          }
      }
    // little-endian limbs: chunk t is the 24 bits at byte 3t
    static_assert(chunk_bits == 24, "byte extraction assumes 24-bit chunks");
    const unsigned char *bytes
      = reinterpret_cast<const unsigned char *>(buffer.data());
    for(int t = 0; t < num_chunks; ++t)
      {
        uint32_t v;
        std::memcpy(&v, bytes + 3 * t, sizeof(v));
        chunks[t * stride] = sign * double(v & 0xFFFFFFu);
      }
  }

  Context::Context(const long bits_A, const long bits_B, const slong k)
      : comb(bits_A + bits_B + FLINT_BIT_COUNT(k) + 3, k),
        num_primes(comb.num_primes),
        num_chunks((std::max(bits_A, bits_B) + chunk_bits - 1) / chunk_bits)
  {
    const long bits = std::min(bits_A, bits_B);
    primes.resize(num_primes);
    inverse_primes.resize(num_primes);
    u.resize(num_primes);
    W.resize(size_t(num_chunks) * num_primes);

    fmpz_t M, M_p;
    fmpz_init_set_ui(M, 1);
    fmpz_init(M_p);
    for(size_t i = 0; i < num_primes; ++i)
      fmpz_mul_ui(M, M, comb.primes[i]);

    // t0 <= (bits - 117) / 24 bounds the error (see blas_residues.hxx).
    // A multiple of 8 makes 24 t0 a multiple of 64: the result is then
    // aligned to limbs as the exact integer, so that its truncation to mpf
    // is at the same bit as for FLINT's CRT. The results then usually
    // agree bit for bit, but not always (see blas_residues.hxx).
    t0 = std::max<long>(0, (bits - 117) / chunk_bits);
    t0 -= t0 % 8;
    const long M_bits = fmpz_bits(M);
    num_high_chunks = (M_bits + chunk_bits - 1) / chunk_bits - t0;
    M_chunks.resize(num_primes * num_high_chunks);
    // S_high < num_primes * 2^53 * 2^(24 num_high_chunks)
    num_limbs = (chunk_bits * num_high_chunks + 64 + 8) / GMP_NUMB_BITS + 2;

    mpz_t z;
    mpz_init(z);
    for(size_t i = 0; i < num_primes; ++i)
      {
        const mp_limb_t p = comb.primes[i];
        ASSERT(std::ldexp(double(p) * num_chunks, chunk_bits)
                   < std::ldexp(1.0, 53)
                 && std::ldexp(double(p) * num_primes, chunk_bits)
                      < std::ldexp(1.0, 53),
               "chunked residues/CRT would not be exact: p=", p,
               " num_chunks=", num_chunks, " num_primes=", num_primes);
        primes[i] = p;
        inverse_primes[i] = 1.0 / double(p);
        mp_limb_t power = 1;
        const mp_limb_t base = (mp_limb_t(1) << chunk_bits) % p;
        for(int t = 0; t < num_chunks; ++t)
          {
            W[t + size_t(num_chunks) * i] = double(power);
            power = (power * base) % p;
          }
        fmpz_divexact_ui(M_p, M, p);
        u[i] = double(n_invmod(fmpz_fdiv_ui(M_p, p), p));
        fmpz_get_mpz(z, M_p);
        for(int t = 0; t < num_high_chunks; ++t)
          {
            // bits [24 (t0+t), 24 (t0+t+1)) of M/p
            mp_limb_t chunk = 0;
            for(int b = 0; b < chunk_bits; ++b)
              chunk |= mp_limb_t(mpz_tstbit(z, chunk_bits * (t0 + t) + b))
                       << b;
            M_chunks[i + num_primes * t] = double(chunk);
          }
      }
    qM.assign((num_primes + 1) * num_limbs, 0);
    for(size_t q = 0; q <= num_primes; ++q)
      {
        fmpz_get_mpz(z, M);
        mpz_mul_ui(z, z, q);
        mpz_fdiv_q_2exp(z, z, chunk_bits * t0);
        ASSERT(mpz_size(z) <= size_t(num_limbs));
        for(size_t l = 0; l < mpz_size(z); ++l)
          qM[q * num_limbs + l] = mpz_getlimbn(z, l);
      }
    mpz_clear(z);
    fmpz_clear(M_p);
    fmpz_clear(M);
  }

  bool is_exact(const long bits_A, const long bits_B, const slong k)
  {
    const Fmpz_Comb comb(bits_A + bits_B + FLINT_BIT_COUNT(k) + 3, k);
    const int num_chunks
      = (std::max(bits_A, bits_B) + chunk_bits - 1) / chunk_bits;
    for(const mp_limb_t p : comb.primes)
      if(!(std::ldexp(double(p) * num_chunks, chunk_bits)
             < std::ldexp(1.0, 53)
           && std::ldexp(double(p) * comb.num_primes, chunk_bits)
                < std::ldexp(1.0, 53)))
        return false;
    return true;
  }

  const Context &
  get_context(const long bits_A, const long bits_B, const slong k)
  {
    static std::map<std::tuple<long, long, slong>, std::unique_ptr<Context>>
      cache;
    auto &context = cache[{bits_A, bits_B, k}];
    if(!context)
      context = std::make_unique<Context>(bits_A, bits_B, k);
    return *context;
  }

  Crt_Workspace::Crt_Workspace(const Context &context)
      : crt_input(context.num_primes), acc(context.num_limbs)
  { fmpz_init(z); }
  Crt_Workspace::~Crt_Workspace()
  { fmpz_clear(z); }

  void blas_crt_high_chunks(const Context &context, const int32_t *residues,
                            const size_t residue_stride, const size_t n,
                            Crt_Workspace &work)
  {
    const size_t num_primes = context.num_primes;
    const int H = context.num_high_chunks;
    work.Y.resize(n * num_primes);
    work.S.resize(n * H);
    for(size_t p = 0; p < num_primes; ++p)
      {
        const int32_t *__restrict y = residues + p * residue_stride;
        double *__restrict yd = work.Y.data() + p * n;
        for(size_t e = 0; e < n; ++e)
          yd[e] = y[e];
      }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, H, num_primes,
                1.0, work.Y.data(), n, context.M_chunks.data(), num_primes,
                0.0, work.S.data(), n);
  }

  namespace
  {
    // c = sum_t S_t 2^(24 t) - q_M, two's complement in acc[0..num_limbs).
    // Returns the sign and sets n = number of significant limbs of |c|.
    int assemble(const Context &context, const double *S, const size_t stride,
                 const mp_limb_t *q_M, mp_limb_t *acc, int &n)
    {
      const int L = context.num_limbs;
      for(int l = 0; l < L; ++l)
        acc[l] = 0;
      for(int t = 0; t < context.num_high_chunks; ++t)
        {
          const mp_limb_t v = mp_limb_t(S[t * stride]);
          const int o = chunk_bits * t;
          int li = o / GMP_NUMB_BITS;
          const int sh = o % GMP_NUMB_BITS;
          mp_limb_t lo = v << sh;
          mp_limb_t hi = sh == 0 ? 0 : v >> (GMP_NUMB_BITS - sh);
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
      mpn_sub_n(acc, acc, q_M, L);
      int sign = 1;
      if(acc[L - 1] >> (GMP_NUMB_BITS - 1))
        {
          sign = -1;
          mpn_neg(acc, acc, L);
        }
      n = L;
      while(n > 0 && acc[n - 1] == 0)
        --n;
      return sign;
    }
  }

  bool blas_crt_to_mpf(const Context &context, const double *S,
                       const size_t stride, const double quotient, long shift,
                       Crt_Workspace &work, mpf_t value)
  {
    const size_t q = std::lround(quotient);
    int size;
    const int sign
      = assemble(context, S, stride, context.qM.data() + q * context.num_limbs,
                 work.acc.data(), size);
    if(size == 0)
      return false;
    __mpf_struct view;
    view._mp_prec = size;
    view._mp_size = sign * size;
    view._mp_exp = size;
    view._mp_d = work.acc.data();
    shift += chunk_bits * context.t0;
    if(shift >= 0)
      mpf_mul_2exp(value, &view, shift);
    else
      mpf_div_2exp(value, &view, -shift);
    return true;
  }

  bool flint_crt_to_mpf(const Context &context, const int32_t *residues,
                        const size_t residue_stride, const long shift,
                        Crt_Workspace &work, mpf_t value)
  {
    // fmpz_multi_CRT_ui writes only to the scratch space comb.comb_temp
    Fmpz_Comb &comb = const_cast<Fmpz_Comb &>(context.comb);
    for(size_t p = 0; p < context.num_primes; ++p)
      work.crt_input[p] = residues[p * residue_stride];
    fmpz_multi_CRT_ui(work.z, work.crt_input.data(), comb.comb, comb.comb_temp,
                      1);
    if(fmpz_is_zero(work.z))
      return false;
    fmpz_get_mpf(value, work.z);
    if(shift >= 0)
      mpf_mul_2exp(value, value, shift);
    else
      mpf_div_2exp(value, value, -shift);
    return true;
  }
}
