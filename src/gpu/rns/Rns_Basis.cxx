// GPU backend for SDPB: RNS basis (host side).  See Rns_Basis.hxx.
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// sdpvibe/rns.py (prime choice, CRT digit tables); portions adapted from SDPviBe,
// Copyright (c) 2026 Antoine Tilloy, MIT License (docs/GPU_backend_credits.md).

#include "Rns_Basis.hxx"

#include <gmp.h>

#include <cmath>
#include <stdexcept>

namespace sdpb_gpu
{
  namespace
  {
    constexpr int prime_max = 16381; // largest prime below 2^14
    constexpr int prime_min = 8192;  // keep primes in (2^13, 2^14)

    std::vector<int32_t> primes_descending()
    {
      std::vector<char> composite(prime_max + 1, 0);
      for(int p = 2; p * p <= prime_max; ++p)
        if(!composite[p])
          for(int q = p * p; q <= prime_max; q += p)
            composite[q] = 1;
      std::vector<int32_t> out;
      for(int p = prime_max; p > prime_min; --p)
        if(!composite[p])
          out.push_back(p);
      return out;
    }

    int centred(long v, long m)
    {
      v %= m;
      if(v < 0)
        v += m;
      return v > m / 2 ? int(v - m) : int(v);
    }

    // balanced base-256 (or base-`base`) digits of x (x is consumed); throws
    // if more are needed
    std::vector<int> balanced_digits(mpz_t x, int count, long base = 256)
    {
      std::vector<int> d(count);
      mpz_t q;
      mpz_init(q);
      for(int j = 0; j < count; ++j)
        {
          // r = x mod base in [0, base), then shift to [-base/2, base/2)
          unsigned long r = mpz_fdiv_q_ui(q, x, base);
          long dr = long(r);
          if(dr >= base / 2)
            {
              dr -= base;
              mpz_add_ui(q, q, 1);
            }
          d[j] = int(dr);
          mpz_swap(x, q);
        }
      if(mpz_sgn(x) != 0)
        throw std::runtime_error("balanced_digits: value needs more digits");
      mpz_clear(q);
      return d;
    }
  }

  Rns_Basis Rns_Basis::for_bound(int bound_bits)
  {
    static const std::vector<int32_t> all = primes_descending();
    Rns_Basis b;
    mpz_t M;
    mpz_init_set_ui(M, 1);
    for(int32_t p : all)
      {
        if(int(mpz_sizeinbase(M, 2)) > bound_bits + 1)
          break;
        b.primes.push_back(p);
        mpz_mul_ui(M, M, p);
      }
    if(int(mpz_sizeinbase(M, 2)) <= bound_bits + 1)
      throw std::runtime_error("Rns_Basis: not enough 14-bit primes");
    b.M_bits = int(mpz_sizeinbase(M, 2));
    b.num_padded = (int(b.primes.size()) + 15) / 16 * 16;
    // CRT sums: |sum_i y_i M/m_i| < num_primes * M / 2, plus the k*M term
    // and a carry digit
    int sum_bits = b.M_bits + int(std::ceil(std::log2(double(b.primes.size())))) + 2;
    b.crt_digits = ((sum_bits + 7) / 8 + 2 + 3) / 4 * 4;
    mpz_clear(M);
    return b;
  }

  Rns_Basis Rns_Basis::for_products(int digits, int64_t rows)
  {
    // |a| < 2^(8 digits - 1)  =>  |sum of rows products| < rows 2^(16 digits - 2);
    // M > 4 |sum| = rows 2^(16 digits)
    int row_bits = int(std::ceil(std::log2(double(rows > 1 ? rows : 2))));
    return for_bound(16 * digits + row_bits);
  }

  void Rns_Basis::residue_tables(int digits, std::vector<int8_t> &T_hi,
                                 std::vector<int8_t> &T_lo) const
  {
    T_hi.assign(size_t(num_padded) * digits, 0);
    T_lo.assign(size_t(num_padded) * digits, 0);
    for(size_t i = 0; i < primes.size(); ++i)
      {
        long m = primes[i], t = 1;
        for(int j = 0; j < digits; ++j)
          {
            int c = centred(t, m); // |c| <= 8190
            int lo = ((c + 64) % 128 + 128) % 128 - 64;
            int hi = (c - lo) / 128;
            T_hi[i * digits + j] = int8_t(hi);
            T_lo[i * digits + j] = int8_t(lo);
            t = t * 256 % m;
          }
      }
  }

  std::vector<int32_t> Rns_Basis::crt_inverses() const
  {
    std::vector<int32_t> inv(num_padded, 0);
    mpz_t M, Mi, r, mod;
    mpz_inits(M, Mi, r, mod, nullptr);
    mpz_set_ui(M, 1);
    for(int32_t p : primes)
      mpz_mul_ui(M, M, p);
    for(size_t i = 0; i < primes.size(); ++i)
      {
        mpz_divexact_ui(Mi, M, primes[i]);
        mpz_set_ui(mod, primes[i]);
        mpz_mod(r, Mi, mod);
        if(!mpz_invert(r, r, mod))
          throw std::runtime_error("crt_inverses: not invertible");
        inv[i] = int32_t(mpz_get_ui(r));
      }
    mpz_clears(M, Mi, r, mod, nullptr);
    return inv;
  }

  void Rns_Basis::crt_tables(std::vector<int8_t> &D,
                             std::vector<int32_t> &M_digits) const
  {
    const int E = crt_digits;
    D.assign(size_t(E) * num_padded, 0);
    mpz_t M, Mi;
    mpz_inits(M, Mi, nullptr);
    mpz_set_ui(M, 1);
    for(int32_t p : primes)
      mpz_mul_ui(M, M, p);
    for(size_t i = 0; i < primes.size(); ++i)
      {
        mpz_divexact_ui(Mi, M, primes[i]);
        auto d = balanced_digits(Mi, E);
        for(int j = 0; j < E; ++j)
          D[size_t(j) * num_padded + i] = int8_t(d[j]);
      }
    auto md = balanced_digits(M, E);
    M_digits.assign(md.begin(), md.end());
    mpz_clears(M, Mi, nullptr);
  }

  int Rns_Basis::crt_digits_128() const
  {
    int sum_bits = M_bits + int(std::ceil(std::log2(double(primes.size())))) + 2;
    // + 1 digit for the high halves' shift, + 2 slack
    return (sum_bits + 6) / 7 + 3;
  }

  void Rns_Basis::crt_tables_128(std::vector<int8_t> &D2,
                                 std::vector<int32_t> &M_digits) const
  {
    const int E2 = crt_digits_128(), rp = num_padded;
    D2.assign(size_t(E2) * 2 * rp, 0);
    mpz_t M, Mi;
    mpz_inits(M, Mi, nullptr);
    mpz_set_ui(M, 1);
    for(int32_t p : primes)
      mpz_mul_ui(M, M, p);
    for(size_t i = 0; i < primes.size(); ++i)
      {
        mpz_divexact_ui(Mi, M, primes[i]);
        auto d = balanced_digits(Mi, E2, 128);
        for(int j = 0; j < E2; ++j)
          {
            D2[size_t(j) * 2 * rp + i] = int8_t(d[j]);
            if(j + 1 < E2)
              D2[size_t(j + 1) * 2 * rp + rp + i] = int8_t(d[j]);
            else if(d[j] != 0)
              throw std::runtime_error("crt_tables_128: top digit");
          }
      }
    auto md = balanced_digits(M, E2, 128);
    M_digits.assign(md.begin(), md.end());
    mpz_clears(M, Mi, nullptr);
  }
}
