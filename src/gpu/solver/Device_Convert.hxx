// Host conversions shared by the PMP and SDPA solver drivers.
#pragma once
#include "Host_Float.hxx"
#include <El.hpp>
#include <algorithm>
#include <cstdlib>
#include <thread>
#include <vector>

namespace sdpb_gpu
{
  constexpr int max_limbs = 79;
  constexpr int64_t MPF_ZERO_EXP = -(int64_t(1) << 40);

  // v -> device float: exactly NB = 32 L - 2 bits, rounded half up on the magnitude
  inline void to_device(const El::BigFloat &v, int L, int32_t *m, int64_t plane,
                        int64_t idx, int64_t *e)
  {
    const __mpf_struct *f = v.gmp_float.get_mpf_t();
    const int size = f->_mp_size;
    if(size == 0)
      {
        for(int t = 0; t < L; ++t)
          m[t * plane + idx] = 0;
        e[idx] = MPF_ZERO_EXP;
        return;
      }
    const int n = size < 0 ? -size : size;
    const mp_limb_t *d = f->_mp_d;
    const int NB = 32 * L - 2;
    const int64_t bits = 64 * int64_t(n) - __builtin_clzll(d[n - 1]);
    int64_t sh = bits - NB;
    uint32_t r[max_limbs + 1];
    auto limb32 = [&](int64_t i) -> uint32_t {
      if(i < 0 || i >= 2 * int64_t(n))
        return 0;
      const mp_limb_t w = d[i / 2];
      return (i & 1) ? uint32_t(w >> 32) : uint32_t(w);
    };
    auto extract = [&](int64_t s) {
      for(int t = 0; t <= L; ++t)
        {
          const int64_t bit = 32 * int64_t(t) + s;
          const int64_t q = bit >= 0 ? bit / 32 : -((-bit + 31) / 32);
          const int rb = int(bit - q * 32);
          r[t] = uint32_t((((uint64_t)limb32(q + 1) << 32) | limb32(q)) >> rb);
        }
      if(s >= 1)
        {
          const int64_t gb = s - 1;
          uint64_t c = (limb32(gb / 32) >> (gb % 32)) & 1u;
          for(int t = 0; t <= L && c; ++t)
            {
              const uint64_t w = uint64_t(r[t]) + c;
              r[t] = uint32_t(w);
              c = w >> 32;
            }
        }
    };
    extract(sh);
    if(r[L] != 0 || (r[L - 1] >> (NB - 32 * (L - 1))) != 0)
      {
        ++sh;
        extract(sh);
      }
    if(size < 0)
      {
        uint64_t c = 1;
        for(int t = 0; t < L; ++t)
          {
            const uint64_t w = uint64_t(~r[t]) + c;
            r[t] = uint32_t(w);
            c = w >> 32;
          }
      }
    for(int t = 0; t < L; ++t)
      m[t * plane + idx] = int32_t(r[t]);
    // the integer of `bits` bits has exponent 64 (exp - n)
    e[idx] = 64 * (int64_t(f->_mp_exp) - n) + sh;
  }

  inline void one_device(int L, int32_t *m, int64_t plane, int64_t idx,
                         int64_t *e)
  {
    const int NB = 32 * L - 2;
    for(int t = 0; t < L; ++t)
      m[t * plane + idx] = t == L - 1 ? int32_t(1u << (NB - 1 - 32 * (L - 1))) : 0;
    e[idx] = -(NB - 1);
  }

  inline void zero_device(int L, int32_t *m, int64_t plane, int64_t idx,
                          int64_t *e)
  {
    for(int t = 0; t < L; ++t)
      m[t * plane + idx] = 0;
    e[idx] = MPF_ZERO_EXP;
  }

  // limbs (two's complement) * 2^e -> BigFloat, exact when it fits the precision
  inline El::BigFloat from_limbs(const uint32_t *limbs, int k, int64_t stride, int64_t e)
  {
    std::vector<uint32_t> a(k);
    for(int t = 0; t < k; ++t)
      a[t] = limbs[t * stride];
    const bool neg = int32_t(a[k - 1]) < 0;
    if(neg)
      {
        uint64_t c = 1;
        for(int t = 0; t < k; ++t)
          {
            const uint64_t w = uint64_t(~a[t]) + c;
            a[t] = uint32_t(w);
            c = w >> 32;
          }
      }
    mpz_t z;
    mpz_init(z);
    mpz_import(z, k, -1, 4, 0, 0, a.data());
    El::BigFloat r;
    mpf_set_z(r.gmp_float.get_mpf_t(), z);
    mpz_clear(z);
    if(e >= 0)
      mpf_mul_2exp(r.gmp_float.get_mpf_t(), r.gmp_float.get_mpf_t(), e);
    else
      mpf_div_2exp(r.gmp_float.get_mpf_t(), r.gmp_float.get_mpf_t(), -e);
    if(neg)
      mpf_neg(r.gmp_float.get_mpf_t(), r.gmp_float.get_mpf_t());
    return r;
  }

  inline El::BigFloat from_host(const Host_Float &h)
  {
    if(h.m.empty() || h.e == MPF_ZERO_EXP)
      return El::BigFloat(0);
    return from_limbs(h.m.data(), int(h.m.size()), 1, h.e);
  }

  inline Host_Float to_host(const El::BigFloat &v, int L)
  {
    Host_Float h;
    std::vector<int32_t> m(L);
    int64_t e;
    to_device(v, L, m.data(), 1, 0, &e);
    h.m.assign(reinterpret_cast<uint32_t *>(m.data()),
               reinterpret_cast<uint32_t *>(m.data()) + L);
    h.e = e;
    return h;
  }

  // fn(i) for 0 <= i < n on OMP_NUM_THREADS threads. Not OpenMP: Elemental's
  // headers change with _OPENMP, and they must read alike in all of SDPB.
  template <class Fn> void parallel_for(int64_t n, const Fn &fn)
  {
    const char *s = getenv("OMP_NUM_THREADS");
    int nt = s ? atoi(s) : int(std::thread::hardware_concurrency());
    nt = int(std::max<int64_t>(
      1, std::min<int64_t>(nt, std::max<int64_t>(1, n / 1024))));
    if(nt == 1)
      {
        for(int64_t i = 0; i < n; ++i)
          fn(i);
        return;
      }
    std::vector<std::thread> pool;
    for(int t = 0; t < nt; ++t)
      pool.emplace_back([&, t]() {
        for(int64_t i = n * t / nt; i < n * (t + 1) / nt; ++i)
          fn(i);
      });
    for(auto &th : pool)
      th.join();
  }

  inline long env_long(const char *name, long dflt)
  {
    const char *s = getenv(name);
    return s ? atol(s) : dflt;
  }
  inline double env_double(const char *name, double dflt)
  {
    const char *s = getenv(name);
    return s ? atof(s) : dflt;
  }
}
