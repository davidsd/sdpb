// GPU backend for SDPB: see Mpf_Ops.hxx.  Based on SDPviBe by Antoine Tilloy
// (https://github.com/AntoineTilloy/SDPviBe), the prototype of this backend:
// sdpvibe/mpf.py and sdpvibe/solver.py.  Portions adapted from SDPviBe,
// Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).

#include "Mpf_Ops.hxx"
#include "../linalg/leaf_fused.cuh"
#include "../rns/Limb_Sizes.hxx"
#include "../rns/rns_kernels.cuh"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

#define OPS_CHECK(x)                                                          \
  do                                                                          \
    {                                                                         \
      cudaError_t e_ = (x);                                                   \
      if(e_ != cudaSuccess)                                                   \
        {                                                                     \
          std::ostringstream ss_;                                             \
          ss_ << __FILE__ << ":" << __LINE__ << ": "                          \
              << cudaGetErrorString(e_);                                      \
          throw std::runtime_error(ss_.str());                                \
        }                                                                     \
  } while(0)

namespace sdpb_gpu
{
  namespace
  {
    typedef uint32_t u32;
    typedef unsigned long long u64;
    typedef long long i64;

    constexpr int TH = 256;
    int64_t nblk(int64_t n) { return (n + TH - 1) / TH; }

    void check_L(int L)
    {
      if(!is_compiled_L(L))
        throw std::invalid_argument(
          "Mpf_Ops: no kernels for this L (rns/Limb_Sizes.hxx)");
    }

    __device__ __forceinline__ i64 entry_index()
    {
      return blockIdx.x * (i64)blockDim.x + threadIdx.x;
    }

    template <int N>
    __device__ __forceinline__ void negate_reg(u32 (&a)[N])
    {
      u64 c = 1;
#pragma unroll
      for(int k = 0; k < N; ++k)
        {
          const u64 t = (u64)(~a[k]) + c;
          a[k] = u32(t);
          c = t >> 32;
        }
    }

    template <int N>
    __device__ __forceinline__ int bitlen_reg(const u32 (&a)[N])
    {
      int bl = 0;
#pragma unroll
      for(int k = 0; k < N; ++k)
        if(a[k])
          bl = 32 * k + 32 - __clz(a[k]);
      return bl;
    }

    // a >>= s for a magnitude, rounding half up, 0 <= s; indices static (barrel
    // over limbs, funnel for the bits), so that everything stays in registers
    template <int N>
    __device__ __forceinline__ void shift_right_round(u32 (&a)[N], i64 s)
    {
      if(s <= 0)
        return;
      if(s >= 32 * N + 1)
        {
#pragma unroll
          for(int k = 0; k < N; ++k)
            a[k] = 0;
          return;
        }
      const int gq = int((s - 1) >> 5), gr = int((s - 1) & 31);
      u32 g = 0;
#pragma unroll
      for(int k = 0; k < N; ++k)
        if(k == gq)
          g = (a[k] >> gr) & 1u;
      const int qs = int(s >> 5), rr = int(s & 31);
#pragma unroll
      for(int b = 0; (1 << b) <= N; ++b)
        if(qs & (1 << b))
          {
#pragma unroll
            for(int k = 0; k < N; ++k)
              a[k] = (k + (1 << b) < N) ? a[k + (1 << b)] : 0u;
          }
#pragma unroll
      for(int k = 0; k < N; ++k)
        a[k] = __funnelshift_r(a[k], (k + 1 < N) ? a[k + 1] : 0u, rr);
      u64 c = g;
#pragma unroll
      for(int k = 0; k < N; ++k)
        {
          const u64 v = (u64)a[k] + c;
          a[k] = u32(v);
          c = v >> 32;
        }
    }

    // normalize a magnitude (W limbs) times 2^base to r (L limbs, exactly NB
    // bits) and exponent: returns false for zero.  r has a round-half-up
    // rounding; a rounding carry into bit NB moves the exponent.
    template <int L, int W>
    __device__ __forceinline__ bool
    round_to_float(const u32 (&a)[W], i64 base, u32 (&r)[L], i64 &ex)
    {
      constexpr int NB = 32 * L - 2;
      const int bl = bitlen_reg<W>(a);
      if(bl == 0)
        return false;
      const int sh = bl - NB;
      if(sh <= 0)
        {
          // left shift: fewer bits than a float holds (no rounding); a fits in
          // L limbs here.  Static indices, so that everything stays in registers
          const int t = -sh, qs = t >> 5, rb = t & 31;
#pragma unroll
          for(int k = 0; k < L; ++k)
            r[k] = a[k];
#pragma unroll
          for(int b = 0; (1 << b) < L; ++b)
            if(qs & (1 << b))
              {
#pragma unroll
                for(int k = L - 1; k >= 0; --k)
                  r[k] = (k - (1 << b) >= 0) ? r[k - (1 << b)] : 0u;
              }
          if(rb)
            {
#pragma unroll
              for(int k = L - 1; k >= 0; --k)
                r[k] = (r[k] << rb) | (k > 0 ? (r[k - 1] >> (32 - rb)) : 0u);
            }
          ex = base + sh;
          return true;
        }
      const int gq = (sh - 1) >> 5, gr = (sh - 1) & 31;
      u32 g = 0;
#pragma unroll
      for(int k = 0; k < W; ++k)
        if(k == gq)
          g = (a[k] >> gr) & 1u;
      u32 w[W];
#pragma unroll
      for(int k = 0; k < W; ++k)
        w[k] = a[k];
      const int qs = sh >> 5, rr = sh & 31;
#pragma unroll
      for(int b = 0; (1 << b) < W; ++b)
        if(qs & (1 << b))
          {
#pragma unroll
            for(int k = 0; k < W; ++k)
              w[k] = (k + (1 << b) < W) ? w[k + (1 << b)] : 0u;
          }
      u32 top = 0;
      u64 c = g;
#pragma unroll
      for(int k = 0; k <= L; ++k)
        {
          const u32 digit
            = __funnelshift_r(w[k], (k + 1 < W) ? w[k + 1] : 0u, rr);
          const u64 v = (u64)digit + c;
          if(k < L)
            r[k] = u32(v);
          else
            top = u32(v);
          c = v >> 32;
        }
      ex = base + sh;
      if(top != 0 || (r[L - 1] >> (NB - 32 * (L - 1))) != 0)
        {
          // 2^NB exactly: one more bit
#pragma unroll
          for(int k = 0; k < L; ++k)
            r[k] = 0;
          r[L - 1] = 1u << (NB - 1 - 32 * (L - 1));
          ex = base + sh + 1;
        }
      return true;
    }

    template <int L>
    __device__ __forceinline__ void read_mag(const int32_t *m, int64_t plane,
                                             int64_t idx, u32 (&a)[L], bool &neg)
    {
#pragma unroll
      for(int t = 0; t < L; ++t)
        a[t] = u32(m[(int64_t)t * plane + idx]);
      neg = int32_t(a[L - 1]) < 0;
      if(neg)
        negate_reg<L>(a);
    }

    template <int L>
    __device__ __forceinline__ void write_float(int32_t *m, int64_t *e,
                                                int64_t plane, int64_t idx,
                                                u32 (&r)[L], i64 ex, bool neg)
    {
      if(neg)
        negate_reg<L>(r);
#pragma unroll
      for(int t = 0; t < L; ++t)
        m[(int64_t)t * plane + idx] = int32_t(r[t]);
      e[idx] = ex;
    }

    template <int L>
    __device__ __forceinline__ void write_zero(int32_t *m, int64_t *e,
                                               int64_t plane, int64_t idx)
    {
#pragma unroll
      for(int t = 0; t < L; ++t)
        m[(int64_t)t * plane + idx] = 0;
      e[idx] = MPF_EZ;
    }

    // (r, ex, neg) = a b, rounded once; magnitudes with exactly NB bits
    template <int L>
    __device__ __forceinline__ void
    mul_float(const u32 (&a)[L], i64 ea, const u32 (&b)[L], i64 eb, u32 (&r)[L],
              i64 &ex)
    {
      // Columns below L - 3 only feed a carry: magnitudes have NB bits, so the
      // product has at least 2 NB - 1 bits, and those columns move it by far
      // less than the rounding position (64 L - 5 - NB = 32 L - 3 bits above bit 0)
      constexpr int KMIN = L - 3, W = 2 * L - KMIN;
      u32 p[W];
      u32 c0 = 0, c1 = 0, c2 = 0;
#pragma unroll
      for(int k = KMIN; k < 2 * L - 1; ++k)
        {
#pragma unroll
          for(int i = (k - (L - 1) > 0 ? k - (L - 1) : 0);
              i <= (k < L - 1 ? k : L - 1); ++i)
            leaf_detail::mac3(c0, c1, c2, a[i], b[k - i]);
          p[k - KMIN] = c0;
          c0 = c1;
          c1 = c2;
          c2 = 0;
        }
      p[W - 1] = c0;
      round_to_float<L, W>(p, ea + eb + 32 * KMIN, r, ex);
    }

    // (r, ex, neg) = (a, ea, na) + (b, eb, nb) for floats
    template <int L>
    __device__ __forceinline__ void
    add_float(const u32 (&a)[L], i64 ea, bool na, const u32 (&b)[L], i64 eb,
              bool nb, u32 (&r)[L], i64 &ex, bool &neg)
    {
      if(ea == MPF_EZ && eb == MPF_EZ)
        {
          ex = MPF_EZ;
          neg = false;
#pragma unroll
          for(int k = 0; k < L; ++k)
            r[k] = 0;
          return;
        }
      constexpr int W = L + 2;
      const i64 t = (ea > eb ? ea : eb) - 32;
      u32 wa[W], wb[W];
      auto shifted = [&](const u32(&x)[L], i64 ee, u32(&w)[W]) {
#pragma unroll
        for(int k = 0; k < W; ++k)
          w[k] = 0;
        if(ee == MPF_EZ || ee - t < -i64(32 * (L + 3)))
          return;
        const i64 d = ee - t;
        if(d >= 0)
          {
            // left shift by 0..32
            const int dd = int(d);
#pragma unroll
            for(int k = 0; k < W; ++k)
              {
                const u32 lo = k < L ? x[k] : 0u;
                const u32 lo2 = (k >= 1 && k - 1 < L) ? x[k - 1] : 0u;
                if(dd == 32)
                  w[k] = lo2;
                else
                  w[k] = dd ? ((lo << dd) | (lo2 >> (32 - dd))) : lo;
              }
          }
        else
          {
#pragma unroll
            for(int k = 0; k < L; ++k)
              w[k] = x[k];
            shift_right_round<W>(w, -d);
          }
      };
      shifted(a, ea, wa);
      shifted(b, eb, wb);
      if(na)
        negate_reg<W>(wa);
      if(nb)
        negate_reg<W>(wb);
      u64 c = 0;
#pragma unroll
      for(int k = 0; k < W; ++k)
        {
          const u64 v = (u64)wa[k] + wb[k] + c;
          wa[k] = u32(v);
          c = v >> 32;
        }
      neg = int32_t(wa[W - 1]) < 0;
      if(neg)
        negate_reg<W>(wa);
      if(!round_to_float<L, W>(wa, t, r, ex))
        {
          ex = MPF_EZ;
          neg = false;
#pragma unroll
          for(int k = 0; k < L; ++k)
            r[k] = 0;
        }
    }

    template <int L>
    __global__ void k_mul(const int32_t *__restrict__ Am,
                          const int64_t *__restrict__ Ae,
                          const int32_t *__restrict__ Bm,
                          const int64_t *__restrict__ Be,
                          int32_t *__restrict__ Om, int64_t *__restrict__ Oe,
                          int64_t n, bool bs)
    {
      const int64_t i = entry_index();
      if(i >= n)
        return;
      const int64_t ib = bs ? 0 : i, pb = bs ? 1 : n;
      const i64 ea = Ae[i], eb = Be[ib];
      if(ea == MPF_EZ || eb == MPF_EZ)
        {
          write_zero<L>(Om, Oe, n, i);
          return;
        }
      u32 a[L], b[L], r[L];
      bool na, nb;
      read_mag<L>(Am, n, i, a, na);
      read_mag<L>(Bm, pb, ib, b, nb);
      i64 ex;
      mul_float<L>(a, ea, b, eb, r, ex);
      write_float<L>(Om, Oe, n, i, r, ex, na != nb);
    }

    template <int L>
    __global__ void k_neg(const int32_t *Am, const int64_t *Ae, int32_t *Om,
                          int64_t *Oe, int64_t n)
    {
      const int64_t i = entry_index();
      if(i >= n)
        return;
      u32 a[L];
#pragma unroll
      for(int t = 0; t < L; ++t)
        a[t] = u32(Am[(int64_t)t * n + i]);
      negate_reg<L>(a);
#pragma unroll
      for(int t = 0; t < L; ++t)
        Om[(int64_t)t * n + i] = int32_t(a[t]);
      Oe[i] = Ae[i];
    }

    __global__ void k_shift_exp(int64_t *e, int64_t n, int64_t k)
    {
      const int64_t i = entry_index();
      if(i < n && e[i] != MPF_EZ)
        e[i] += k;
    }

    template <int L>
    __global__ void k_lower_to_upper(int32_t *m, int64_t *e, int64_t n,
                                     int batch)
    {
      const int64_t q = entry_index(), tot = int64_t(batch) * n * n;
      if(q >= tot)
        return;
      const int64_t b = q / (n * n), rem = q % (n * n);
      const int64_t i = rem % n, j = rem / n;
      if(i >= j)
        return;
      const int64_t from = b * n * n + i * n + j;
      e[q] = e[from];
#pragma unroll
      for(int t = 0; t < L; ++t)
        m[(int64_t)t * tot + q] = m[(int64_t)t * tot + from];
    }

    template <int L>
    __global__ void k_mask_blocks(int32_t *m, int64_t *e, int64_t n, int batch,
                                  int h, const int *hj)
    {
      const int64_t q = entry_index(), tot = int64_t(batch) * n * n;
      if(q >= tot)
        return;
      const int64_t b = q / (n * n), rem = q % (n * n);
      const int i = int(rem % n), j = int(rem / n), hb = hj[b];
      if((i % h) < hb && (j % h) < hb)
        return;
      write_zero<L>(m, e, tot, q);
    }

    template <int L>
    __global__ void
    k_add_diag(int32_t *__restrict__ m, int64_t *__restrict__ e, int64_t n,
               int batch, const int32_t *__restrict__ am,
               const int64_t *__restrict__ ae, int h, const int *__restrict__ hj)
    {
      const int64_t q = entry_index();
      if(q >= int64_t(batch) * n)
        return;
      const int64_t b = q / n, i = q % n;
      if(int(i % h) >= hj[b])
        return;
      const int64_t tot = int64_t(batch) * n * n, idx = b * n * n + i * n + i;
      u32 x[L], s[L], r[L];
      bool nx, ns, neg;
      i64 ex;
      read_mag<L>(m, tot, idx, x, nx);
      read_mag<L>(am, 1, 0, s, ns);
      add_float<L>(x, e[idx], nx, s, ae[0], ns, r, ex, neg);
      write_float<L>(m, e, tot, idx, r, ex, neg);
    }

    template <int L>
    __global__ void k_pad_identity(const int32_t *__restrict__ sm,
                                   const int64_t *__restrict__ se, int64_t n,
                                   int32_t *__restrict__ dm,
                                   int64_t *__restrict__ de, int64_t np,
                                   int batch)
    {
      const int64_t q = entry_index(), tot = int64_t(batch) * np * np;
      if(q >= tot)
        return;
      const int64_t b = q / (np * np), rem = q % (np * np);
      const int64_t i = rem % np, j = rem / np;
      if(i < n && j < n)
        {
          const int64_t from = b * n * n + j * n + i,
                        stot = int64_t(batch) * n * n;
          de[q] = se[from];
#pragma unroll
          for(int t = 0; t < L; ++t)
            dm[(int64_t)t * tot + q] = sm[(int64_t)t * stot + from];
        }
      else if(i == j)
        {
          constexpr int NB = 32 * L - 2;
#pragma unroll
          for(int t = 0; t < L; ++t)
            dm[(int64_t)t * tot + q]
              = t == L - 1 ? int32_t(1u << (NB - 1 - 32 * (L - 1))) : 0;
          de[q] = -(NB - 1);
        }
      else
        write_zero<L>(dm, de, tot, q);
    }

    template <int L>
    __global__ void k_diagonal(const int32_t *__restrict__ m,
                               const int64_t *__restrict__ e, int64_t n,
                               int batch, int32_t *__restrict__ dm,
                               int64_t *__restrict__ de)
    {
      const int64_t q = entry_index(), tot = int64_t(batch) * n;
      if(q >= tot)
        return;
      const int64_t b = q / n, i = q % n;
      const int64_t from = b * n * n + i * n + i,
                    ftot = int64_t(batch) * n * n;
      de[q] = e[from];
#pragma unroll
      for(int t = 0; t < L; ++t)
        dm[(int64_t)t * tot + q] = m[(int64_t)t * ftot + from];
    }

    template <int L>
    __global__ void k_kron(const int32_t *__restrict__ vm,
                           const int64_t *__restrict__ ve, int64_t h, int64_t K,
                           int dim, int batch, int32_t *__restrict__ om,
                           int64_t *__restrict__ oe)
    {
      const int64_t R = dim * h, C = dim * K;
      const int64_t q = entry_index(), tot = int64_t(batch) * R * C;
      if(q >= tot)
        return;
      const int64_t b = q / (R * C), rem = q % (R * C);
      const int64_t i = rem % R, j = rem / R;
      const int64_t ci = i / h, cj = j / K;
      if(ci != cj)
        {
          write_zero<L>(om, oe, tot, q);
          return;
        }
      const int64_t from = b * h * K + (j % K) * h + (i % h),
                    ftot = int64_t(batch) * h * K;
      oe[q] = ve[from];
#pragma unroll
      for(int t = 0; t < L; ++t)
        om[(int64_t)t * tot + q] = vm[(int64_t)t * ftot + from];
    }

    template <int L>
    __global__ void k_repeat(const int32_t *__restrict__ am,
                             const int64_t *__restrict__ ae, int64_t ent,
                             int32_t *__restrict__ om, int64_t *__restrict__ oe,
                             int copies)
    {
      const int64_t q = entry_index(), tot = ent * copies;
      if(q >= tot)
        return;
      const int64_t from = q % ent;
      oe[q] = ae[from];
#pragma unroll
      for(int t = 0; t < L; ++t)
        om[(int64_t)t * tot + q] = am[(int64_t)t * ent + from];
    }

    __device__ __forceinline__ void seg_cr(int s, int &c, int &r)
    {
      c = 0;
      while((c + 1) * (c + 2) / 2 <= s)
        ++c;
      r = s - c * (c + 1) / 2;
    }

    template <int L>
    __global__ void k_seg_diag(const int32_t *__restrict__ mm,
                               const int64_t *__restrict__ me, int64_t D,
                               int nb, int64_t K, int nseg,
                               int32_t *__restrict__ vm, int64_t *__restrict__ ve)
    {
      const int64_t P = nseg * K, q = entry_index(), tot = int64_t(nb) * P;
      if(q >= tot)
        return;
      const int64_t b = q / P, rem = q % P;
      const int s = int(rem / K);
      const int64_t k = rem % K;
      int c, r;
      seg_cr(s, c, r);
      const int64_t ftot = int64_t(2 * nb) * D * D;
      const int64_t off = (c * K + k) * D + r * K + k;
      const int64_t i0 = b * D * D + off, i1 = (b + nb) * D * D + off;
      u32 x[L], y[L], o[L];
      bool nx, ny, neg;
      i64 ex;
      read_mag<L>(mm, ftot, i0, x, nx);
      read_mag<L>(mm, ftot, i1, y, ny);
      add_float<L>(x, me[i0], nx, y, me[i1], ny, o, ex, neg);
      write_float<L>(vm, ve, tot, q, o, ex, neg);
    }

    template <int L>
    __global__ void k_scale_segments(const int32_t *__restrict__ vm,
                                     const int64_t *__restrict__ ve, int64_t K,
                                     int64_t h, int nb,
                                     const int32_t *__restrict__ am,
                                     const int64_t *__restrict__ ae, int nseg,
                                     int32_t *__restrict__ om,
                                     int64_t *__restrict__ oe)
    {
      const int64_t ent = K * h, batch = int64_t(2 * nb) * nseg,
                    tot = batch * ent;
      const int64_t q = entry_index();
      if(q >= tot)
        return;
      const int64_t bb = q / ent, rem = q % ent;
      const int64_t z = bb % (2 * nb), s = bb / (2 * nb), qq = z % nb;
      const int64_t k = rem % K;
      const int64_t vidx = z * ent + rem, ftot = int64_t(2 * nb) * ent;
      const int64_t P = nseg * K, aidx = qq * P + s * K + k,
                    atot = int64_t(nb) * P;
      const i64 ev = ve[vidx], ea = ae[aidx];
      if(ev == MPF_EZ || ea == MPF_EZ)
        {
          write_zero<L>(om, oe, tot, q);
          return;
        }
      u32 x[L], y[L], r[L];
      bool nx, ny;
      i64 ex;
      read_mag<L>(vm, ftot, vidx, x, nx);
      read_mag<L>(am, atot, aidx, y, ny);
      mul_float<L>(x, ev, y, ea, r, ex);
      write_float<L>(om, oe, tot, q, r, ex, nx != ny);
    }

    template <int L>
    __global__ void k_cmws(const int32_t *__restrict__ cm,
                           const int64_t *__restrict__ ce, int64_t h, int dim,
                           int nbp, int32_t *__restrict__ om,
                           int64_t *__restrict__ oe)
    {
      const int64_t n = dim * h, q = entry_index(), tot = int64_t(nbp) * n * n;
      if(q >= tot)
        return;
      const int nseg = dim * (dim + 1) / 2;
      const int64_t z = q / (n * n), rem = q % (n * n);
      const int64_t i0 = rem % n, j0 = rem / n;
      const int64_t ri = i0 / h, cj = j0 / h, i = i0 % h, j = j0 % h;
      // the block (ri, cj): for ri <= cj the segment (c = cj, r = ri) at (i, j),
      // else its transpose: the segment (c = ri, r = cj) at (j, i)
      const int64_t c = ri <= cj ? cj : ri, r = ri <= cj ? ri : cj;
      const int64_t s = c * (c + 1) / 2 + r;
      const int64_t a = ri <= cj ? i : j, bcol = ri <= cj ? j : i;
      const int64_t ftot = int64_t(nbp) * nseg * h * h;
      const int64_t from = (z + nbp * s) * h * h + bcol * h + a;
      const i64 ex = ce[from];
      oe[q] = (ex != MPF_EZ && c != r) ? ex - 1 : ex;
#pragma unroll
      for(int t = 0; t < L; ++t)
        om[(int64_t)t * tot + q] = cm[(int64_t)t * ftot + from];
    }

    template <int L>
    __global__ void
    k_schur(const int32_t *__restrict__ xm, const int64_t *__restrict__ xe,
            const int32_t *__restrict__ ym, const int64_t *__restrict__ ye,
            int64_t K, int dim, int nb, const int *__restrict__ Kj,
            int32_t *__restrict__ sm, int64_t *__restrict__ se)
    {
      const int nseg = dim * (dim + 1) / 2;
      const int64_t P = nseg * K, tot = int64_t(nb) * P * P;
      const int64_t q = entry_index();
      if(q >= tot)
        return;
      const int64_t b = q / (P * P), rem = q % (P * P);
      const int64_t row = rem % P, col = rem / P;
      if(row < col)
        return;
      const int64_t mirror = b * P * P + row * P + col;
      const int64_t k1 = row % K, k2 = col % K;
      if(k1 >= Kj[b] || k2 >= Kj[b])
        {
          if(row == col)
            {
              constexpr int NB = 32 * L - 2;
#pragma unroll
              for(int t = 0; t < L; ++t)
                sm[(int64_t)t * tot + q]
                  = t == L - 1 ? int32_t(1u << (NB - 1 - 32 * (L - 1))) : 0;
              se[q] = -(NB - 1);
            }
          else
            {
              write_zero<L>(sm, se, tot, q);
              write_zero<L>(sm, se, tot, mirror);
            }
          return;
        }
      int c0, r0, c1, r1;
      seg_cr(int(row / K), c0, r0);
      seg_cr(int(col / K), c1, r1);
      const int64_t D = dim * K, ftot = int64_t(2 * nb) * D * D;
      // A(u K + k1, v K + k2) in the column-major block z
      auto at = [&](int z, int u, int v) {
        return z * D * D + (v * K + k2) * D + u * K + k1;
      };
      u32 acc[L];
      i64 ea = MPF_EZ;
      bool na = false;
#pragma unroll
      for(int t = 0; t < L; ++t)
        acc[t] = 0;
      const int nterms = dim == 1 ? 1 : 4;
      for(int z = int(b); z <= int(b) + nb; z += nb)
        for(int tn = 0; tn < nterms; ++tn)
          {
            // SDPB's four terms (compute_schur_complement.cxx); for dim = 1 they
            // coincide, one term and no factor 1/4
            const int xu = (tn == 0 || tn == 2) ? c0 : r0;
            const int xv = tn < 2 ? r1 : c1;
            const int yu = (tn == 0 || tn == 2) ? r0 : c0;
            const int yv = tn == 0 ? c1 : (tn == 1 ? c1 : r1);
            const int64_t ix = at(z, xu, xv), iy = at(z, yu, yv);
            const i64 ex = xe[ix], ey = ye[iy];
            if(ex == MPF_EZ || ey == MPF_EZ)
              continue;
            u32 x[L], y[L], p[L], o[L];
            bool nx, ny, neg;
            i64 ep, eo;
            read_mag<L>(xm, ftot, ix, x, nx);
            read_mag<L>(ym, ftot, iy, y, ny);
            mul_float<L>(x, ex, y, ey, p, ep);
            add_float<L>(acc, ea, na, p, ep, nx != ny, o, eo, neg);
#pragma unroll
            for(int t = 0; t < L; ++t)
              acc[t] = o[t];
            ea = eo;
            na = neg;
          }
      if(ea != MPF_EZ && dim != 1)
        ea -= 2;
      if(ea == MPF_EZ)
        {
          write_zero<L>(sm, se, tot, q);
          if(row != col)
            write_zero<L>(sm, se, tot, mirror);
          return;
        }
      if(na)
        negate_reg<L>(acc);
#pragma unroll
      for(int t = 0; t < L; ++t)
        sm[(int64_t)t * tot + q] = int32_t(acc[t]);
      se[q] = ea;
      if(row != col)
        {
#pragma unroll
          for(int t = 0; t < L; ++t)
            sm[(int64_t)t * tot + mirror] = int32_t(acc[t]);
          se[mirror] = ea;
        }
    }

    // A direct product: one thread per output entry.  Terms are aligned to the
    // largest term's exponent T by shifting the matrix entry's magnitude,
    // multiplied (columns below L - 3 only contribute carries, which are
    // dropped, far below the rounding position) and accumulated in two's
    // complement over the columns L - 3 .. 2 L + 1 (2 limbs of headroom).
    template <int L>
    __global__ void
    k_gemv(const int32_t *__restrict__ Am, const int64_t *__restrict__ Ae,
           int64_t aplane, const int32_t *__restrict__ Vm,
           const int64_t *__restrict__ Ve, int64_t vplane, int64_t no,
           int batch, int64_t ab, int64_t ao, int64_t atb, int64_t ati,
           int64_t vb, int64_t vo, int64_t vtb, int64_t vti, int64_t ntb,
           int64_t nti, int64_t ob, int64_t oo, int32_t *__restrict__ Om,
           int64_t *__restrict__ Oe, int64_t oplane)
    {
      const int64_t q = entry_index();
      if(q >= int64_t(batch) * no)
        return;
      const int64_t b = q / no, o = q % no;
      const int64_t a0 = b * ab + o * ao, v0 = b * vb + o * vo;
      const int64_t out = b * ob + o * oo;
      constexpr int NB = 32 * L - 2, KMIN = L - 3, W = 2 * L + 2 - KMIN;
      i64 T = LLONG_MIN;
      for(int64_t tb = 0; tb < ntb; ++tb)
        for(int64_t ti = 0; ti < nti; ++ti)
          {
            const i64 ea = Ae[a0 + tb * atb + ti * ati],
                      ev = Ve[v0 + tb * vtb + ti * vti];
            if(ea != MPF_EZ && ev != MPF_EZ && ea + ev > T)
              T = ea + ev;
          }
      if(T == LLONG_MIN)
        {
          write_zero<L>(Om, Oe, oplane, out);
          return;
        }
      u32 acc[W];
#pragma unroll
      for(int k = 0; k < W; ++k)
        acc[k] = 0;
      for(int64_t tb = 0; tb < ntb; ++tb)
        for(int64_t ti = 0; ti < nti; ++ti)
          {
            const int64_t ia = a0 + tb * atb + ti * ati,
                          iv = v0 + tb * vtb + ti * vti;
            const i64 ea = Ae[ia], ev = Ve[iv];
            if(ea == MPF_EZ || ev == MPF_EZ)
              continue;
            const i64 s = T - (ea + ev);
            if(s > NB)
              continue;
            u32 x[L], y[L];
            bool nx, ny;
            read_mag<L>(Am, aplane, ia, x, nx);
            read_mag<L>(Vm, vplane, iv, y, ny);
            shift_right_round<L>(x, s);
            const bool neg = nx != ny;
            const u32 mask = neg ? 0xffffffffu : 0u;
            u32 cc = neg ? 1u : 0u;
            u32 c0 = 0, c1 = 0, c2 = 0;
#pragma unroll
            for(int k = KMIN; k < 2 * L - 1; ++k)
              {
#pragma unroll
                for(int i = (k - (L - 1) > 0 ? k - (L - 1) : 0);
                    i <= (k < L - 1 ? k : L - 1); ++i)
                  leaf_detail::mac3(c0, c1, c2, x[i], y[k - i]);
                const u64 v = (u64)acc[k - KMIN] + (c0 ^ mask) + cc;
                acc[k - KMIN] = u32(v);
                cc = u32(v >> 32);
                c0 = c1;
                c1 = c2;
                c2 = 0;
              }
            // columns 2 L - 1 .. 2 L + 1: the last carry word, then sign extension
#pragma unroll
            for(int k = 2 * L - 1; k <= 2 * L + 1; ++k)
              {
                const u64 v = (u64)acc[k - KMIN] + (c0 ^ mask) + cc;
                acc[k - KMIN] = u32(v);
                cc = u32(v >> 32);
                c0 = c1;
                c1 = 0;
              }
          }
      const bool neg = int32_t(acc[W - 1]) < 0;
      if(neg)
        negate_reg<W>(acc);
      u32 r[L];
      i64 ex;
      if(round_to_float<L, W>(acc, T + 32 * KMIN, r, ex))
        write_float<L>(Om, Oe, oplane, out, r, ex, neg);
      else
        write_zero<L>(Om, Oe, oplane, out);
    }

    // The same product with several threads per output entry, for shapes with
    // few outputs (Q^-1 v, B^T v, P^T v: N outputs), where one thread per output
    // leaves the streaming multiprocessors nearly idle.  After Zechuan Zheng's
    // gemv in his GPU port of SDPviBe: each of S threads sums
    // every S-th term exactly, in two's complement, and one reduction adds the
    // partial sums and rounds.  All sums are exact and T is the same, so the result
    // is bitwise that of k_gemv.
    template <int L>
    __global__ void
    k_gemv_top(const int64_t *__restrict__ Ae, const int64_t *__restrict__ Ve,
               int64_t no, int batch, int64_t ab, int64_t ao, int64_t atb,
               int64_t ati, int64_t vb, int64_t vo, int64_t vtb, int64_t vti,
               int64_t ntb, int64_t nti, i64 *__restrict__ T)
    {
      const int64_t q = entry_index();
      if(q >= int64_t(batch) * no)
        return;
      const int64_t b = q / no, o = q % no;
      const int64_t a0 = b * ab + o * ao, v0 = b * vb + o * vo;
      i64 top = LLONG_MIN;
      for(int64_t tb = 0; tb < ntb; ++tb)
        for(int64_t ti = 0; ti < nti; ++ti)
          {
            const i64 ea = Ae[a0 + tb * atb + ti * ati],
                      ev = Ve[v0 + tb * vtb + ti * vti];
            if(ea != MPF_EZ && ev != MPF_EZ && ea + ev > top)
              top = ea + ev;
          }
      T[q] = top;
    }

    // lanes share an output and take consecutive terms if by_term (the terms are
    // adjacent in memory), else share a term index and take consecutive outputs
    template <int L>
    __global__ void
    k_gemv_part(const int32_t *__restrict__ Am, const int64_t *__restrict__ Ae,
                int64_t aplane, const int32_t *__restrict__ Vm,
                const int64_t *__restrict__ Ve, int64_t vplane,
                const i64 *__restrict__ T, int64_t no, int batch, int64_t ab,
                int64_t ao, int64_t atb, int64_t ati, int64_t vb, int64_t vo,
                int64_t vtb, int64_t vti, int64_t ntb, int64_t nti, int S,
                int by_term, u32 *__restrict__ Part)
    {
      const int64_t gid = entry_index(), n = int64_t(batch) * no;
      const int64_t q = by_term ? gid / S : gid % n, sl = by_term ? gid % S : gid / n;
      if(q >= n || sl >= S)
        return;
      const int64_t b = q / no, o = q % no;
      const int64_t a0 = b * ab + o * ao, v0 = b * vb + o * vo;
      constexpr int NB = 32 * L - 2, KMIN = L - 3, W = 2 * L + 2 - KMIN;
      const i64 top = T[q];
      u32 acc[W];
#pragma unroll
      for(int k = 0; k < W; ++k)
        acc[k] = 0;
      if(top != LLONG_MIN)
        {
          int64_t first = sl; // the term index in this row of terms, over tb
          for(int64_t tb = 0; tb < ntb; ++tb)
            {
              int64_t ti = first;
              for(; ti < nti; ti += S)
                {
                  const int64_t ia = a0 + tb * atb + ti * ati,
                                iv = v0 + tb * vtb + ti * vti;
                  const i64 ea = Ae[ia], ev = Ve[iv];
                  if(ea == MPF_EZ || ev == MPF_EZ)
                    continue;
                  const i64 s = top - (ea + ev);
                  if(s > NB)
                    continue;
                  u32 x[L], y[L];
                  bool nx, ny;
                  read_mag<L>(Am, aplane, ia, x, nx);
                  read_mag<L>(Vm, vplane, iv, y, ny);
                  shift_right_round<L>(x, s);
                  const bool neg = nx != ny;
                  const u32 mask = neg ? 0xffffffffu : 0u;
                  u32 cc = neg ? 1u : 0u;
                  u32 c0 = 0, c1 = 0, c2 = 0;
#pragma unroll
                  for(int k = KMIN; k < 2 * L - 1; ++k)
                    {
#pragma unroll
                      for(int i = (k - (L - 1) > 0 ? k - (L - 1) : 0);
                          i <= (k < L - 1 ? k : L - 1); ++i)
                        leaf_detail::mac3(c0, c1, c2, x[i], y[k - i]);
                      const u64 v = (u64)acc[k - KMIN] + (c0 ^ mask) + cc;
                      acc[k - KMIN] = u32(v);
                      cc = u32(v >> 32);
                      c0 = c1;
                      c1 = c2;
                      c2 = 0;
                    }
#pragma unroll
                  for(int k = 2 * L - 1; k <= 2 * L + 1; ++k)
                    {
                      const u64 v = (u64)acc[k - KMIN] + (c0 ^ mask) + cc;
                      acc[k - KMIN] = u32(v);
                      cc = u32(v >> 32);
                      c0 = c1;
                      c1 = 0;
                    }
                }
              first = ti - nti;
            }
        }
#pragma unroll
      for(int k = 0; k < W; ++k)
        Part[(sl * W + k) * n + q] = acc[k];
    }

    template <int L>
    __global__ void
    k_gemv_reduce(const u32 *__restrict__ Part, const i64 *__restrict__ T,
                  int64_t no, int batch, int S, int64_t ob, int64_t oo,
                  int32_t *__restrict__ Om, int64_t *__restrict__ Oe,
                  int64_t oplane)
    {
      const int64_t q = entry_index(), n = int64_t(batch) * no;
      if(q >= n)
        return;
      const int64_t out = (q / no) * ob + (q % no) * oo;
      constexpr int KMIN = L - 3, W = 2 * L + 2 - KMIN;
      if(T[q] == LLONG_MIN)
        {
          write_zero<L>(Om, Oe, oplane, out);
          return;
        }
      u32 acc[W];
#pragma unroll
      for(int k = 0; k < W; ++k)
        acc[k] = Part[int64_t(k) * n + q];
      for(int sl = 1; sl < S; ++sl)
        {
          u64 c = 0;
#pragma unroll
          for(int k = 0; k < W; ++k)
            {
              const u64 v = (u64)acc[k] + Part[(int64_t(sl) * W + k) * n + q] + c;
              acc[k] = u32(v);
              c = v >> 32;
            }
        }
      const bool neg = int32_t(acc[W - 1]) < 0;
      if(neg)
        negate_reg<W>(acc);
      u32 r[L];
      i64 ex;
      if(round_to_float<L, W>(acc, T[q] + 32 * KMIN, r, ex))
        write_float<L>(Om, Oe, oplane, out, r, ex, neg);
      else
        write_zero<L>(Om, Oe, oplane, out);
    }

    __global__ void k_max_exp2(const int64_t *e, int64_t n, u64 *out)
    {
      __shared__ i64 v[TH];
      i64 mx = LLONG_MIN;
      for(int64_t i = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; i < n;
          i += int64_t(gridDim.x) * blockDim.x)
        if(e[i] != MPF_EZ)
          mx = max(mx, (i64)e[i]);
      v[threadIdx.x] = mx;
      __syncthreads();
      for(int d = TH / 2; d; d /= 2)
        {
          if(threadIdx.x < d)
            v[threadIdx.x] = max(v[threadIdx.x], v[threadIdx.x + d]);
          __syncthreads();
        }
      if(threadIdx.x == 0)
        atomicMax(out, (u64)v[0] ^ 0x8000000000000000ULL);
    }

    // the top 64 bits of |m| among the entries with e == emax, and an index
    template <int L>
    __global__ void k_max_key(const int32_t *__restrict__ m,
                              const int64_t *__restrict__ e, int64_t n, i64 emax,
                              u64 *__restrict__ key)
    {
      const int64_t i = entry_index();
      if(i >= n || e[i] != emax)
        return;
      u32 a[L];
      bool neg;
      read_mag<L>(m, n, i, a, neg);
      atomicMax(key, ((u64)a[L - 1] << 32) | a[L - 2]);
    }

    // every entry whose top 64 bits are the key goes to the host, which compares
    // all the limbs (entries can share their top 64 bits; a race here would make
    // the reported maximum depend on the run)
    template <int L>
    __global__ void k_find_key(const int32_t *__restrict__ m,
                               const int64_t *__restrict__ e, int64_t n, i64 emax,
                               u64 keyv, int32_t *__restrict__ out,
                               unsigned int *__restrict__ count, unsigned int cap)
    {
      const int64_t i = entry_index();
      if(i >= n || e[i] != emax)
        return;
      u32 a[L];
      bool neg;
      read_mag<L>(m, n, i, a, neg);
      if((((u64)a[L - 1] << 32) | a[L - 2]) != keyv)
        return;
      const unsigned int slot = atomicAdd(count, 1u);
      if(slot >= cap)
        return;
#pragma unroll
      for(int t = 0; t < L; ++t)
        out[slot * L + t] = int32_t(a[t]);
    }

    // limb sums of all entries shifted to the exponent t, block by block
    template <int L>
    __global__ void k_total(const int32_t *__restrict__ m,
                            const int64_t *__restrict__ e, int64_t n, i64 t0,
                            u64 *__restrict__ out)
    {
      constexpr int W = L + 2;
      __shared__ u64 sum[W + 1];
      if(threadIdx.x == 0)
        {
#pragma unroll
          for(int k = 0; k <= W; ++k)
            sum[k] = 0;
        }
      __syncthreads();
      const int64_t i = entry_index();
      if(i < n && e[i] != MPF_EZ && t0 - e[i] <= 32 * L)
        {
          u32 w[W];
          bool neg;
          u32 a[L];
          read_mag<L>(m, n, i, a, neg);
#pragma unroll
          for(int k = 0; k < W; ++k)
            w[k] = k < L ? a[k] : 0u;
          const i64 d = e[i] - t0;
          if(d >= 0)
            {
              // left shift by d <= 32 (t0 = emax - 32)
              const int dd = int(d);
#pragma unroll
              for(int k = W - 1; k >= 0; --k)
                {
                  const u32 lo = w[k], lo2 = k ? w[k - 1] : 0u;
                  w[k] = dd == 0 ? lo
                                 : (dd == 32 ? lo2
                                             : ((lo << dd) | (lo2 >> (32 - dd))));
                }
            }
          else
            shift_right_round<W>(w, -d);
          // a magnitude that rounded to zero is no negative term: -0 would
          // count -2^(32 W) in the sign limb
          u32 any = 0;
#pragma unroll
          for(int k = 0; k < W; ++k)
            any |= w[k];
          if(neg && any)
            negate_reg<W>(w);
#pragma unroll
          for(int k = 0; k < W; ++k)
            atomicAdd(&sum[k], (u64)w[k]);
          // the sign limb
          if(neg && any)
            atomicAdd(&sum[W], 1ULL);
        }
      __syncthreads();
      if(threadIdx.x == 0)
        {
#pragma unroll
          for(int k = 0; k <= W; ++k)
            atomicAdd(&out[k], (u64)sum[k]);
        }
    }

    __global__ void k_block_top(const int64_t *e, int64_t per, int nblocks,
                                u64 *out)
    {
      // batch b = blockIdx.y
      const int64_t b = blockIdx.y;
      __shared__ i64 v[TH];
      i64 mx = LLONG_MIN;
      for(int64_t i = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; i < per;
          i += int64_t(gridDim.x) * blockDim.x)
        {
          const int64_t ev = e[b * per + i];
          if(ev != MPF_EZ)
            mx = max(mx, (i64)ev);
        }
      v[threadIdx.x] = mx;
      __syncthreads();
      for(int d = TH / 2; d; d /= 2)
        {
          if(threadIdx.x < d)
            v[threadIdx.x] = max(v[threadIdx.x], v[threadIdx.x + d]);
          __syncthreads();
        }
      if(threadIdx.x == 0)
        atomicMax(out + b, (u64)v[0] ^ 0x8000000000000000ULL);
    }

    template <int L>
    __global__ void k_to_double(const int32_t *__restrict__ m,
                                const int64_t *__restrict__ e, int64_t per,
                                int batch, const int64_t *__restrict__ top,
                                double *__restrict__ out)
    {
      const int64_t q = entry_index(), tot = int64_t(batch) * per;
      if(q >= tot)
        return;
      const i64 ee = e[q];
      if(ee == MPF_EZ)
        {
          out[q] = 0;
          return;
        }
      const int64_t b = q / per;
      const double hi = double(m[(int64_t)(L - 1) * tot + q]),
                   lo = double(u32(m[(int64_t)(L - 2) * tot + q]));
      i64 sc = ee + 32 * (L - 2) - top[b];
      sc = sc < -1100 ? -1100 : (sc > 1100 ? 1100 : sc);
      out[q] = ldexp(hi * 4294967296.0 + lo, int(sc));
    }
  }

// kernel<L> for the run-time L, among the sizes in rns/Limb_Sizes.hxx
#define OPS_L_CASE(LT, kernel, ...)                                           \
  case LT:                                                                    \
    kernel<LT> __VA_ARGS__;                                                   \
    break;
#define OPS_BY_L(L, kernel, ...)                                              \
  switch(L)                                                                   \
    {                                                                         \
      SDPB_GPU_EACH_L(OPS_L_CASE, kernel, __VA_ARGS__)                        \
    default:                                                                  \
      check_L(L);                                                             \
    }

  void dev_mul(const Mpf_Batch &A, const Mpf_Batch &B, Mpf_Batch &O,
               cudaStream_t s)
  {
    check_L(A.L);
    const int64_t n = A.entries();
    const bool bs = B.entries() == 1;
    if(!bs && B.entries() != n)
      throw std::invalid_argument("dev_mul: shapes");
    OPS_BY_L(A.L, k_mul, <<<nblk(n), TH, 0, s>>>(A.m, A.e, B.m, B.e, O.m, O.e, n, bs));
    OPS_CHECK(cudaGetLastError());
  }

  void dev_neg(const Mpf_Batch &A, Mpf_Batch &O, cudaStream_t s)
  {
    const int64_t n = A.entries();
    OPS_BY_L(A.L, k_neg, <<<nblk(n), TH, 0, s>>>(A.m, A.e, O.m, O.e, n));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_shift_exp(Mpf_Batch &A, int64_t k, cudaStream_t s)
  {
    k_shift_exp<<<nblk(A.entries()), TH, 0, s>>>(A.e, A.entries(), k);
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_lower_to_upper(Mpf_Batch &M, cudaStream_t s)
  {
    OPS_BY_L(M.L, k_lower_to_upper, <<<nblk(M.entries()), TH, 0, s>>>(M.m, M.e, M.rows,
                                                           M.batch));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_mask_blocks(Mpf_Batch &M, int h, const int *d_hj, cudaStream_t s)
  {
    OPS_BY_L(M.L, k_mask_blocks, <<<nblk(M.entries()), TH, 0, s>>>(M.m, M.e, M.rows, M.batch,
                                                        h, d_hj));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_add_diag(Mpf_Batch &M, const Mpf_Batch &a, int h, const int *d_hj,
                    cudaStream_t s)
  {
    OPS_BY_L(M.L, k_add_diag, <<<nblk(int64_t(M.batch) * M.rows), TH, 0, s>>>(
      M.m, M.e, M.rows, M.batch, a.m, a.e, h, d_hj));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_pad_identity(const Mpf_Batch &src, Mpf_Batch &dst, cudaStream_t s)
  {
    OPS_BY_L(src.L, k_pad_identity, <<<nblk(dst.entries()), TH, 0, s>>>(
      src.m, src.e, src.rows, dst.m, dst.e, dst.rows, src.batch));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_diagonal(const Mpf_Batch &M, Mpf_Batch &D, cudaStream_t s)
  {
    OPS_BY_L(M.L, k_diagonal, <<<nblk(D.entries()), TH, 0, s>>>(M.m, M.e, M.rows, M.batch,
                                                     D.m, D.e));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_kron_identity(const Mpf_Batch &V, Mpf_Batch &Vx, int dim,
                         cudaStream_t s)
  {
    OPS_BY_L(V.L, k_kron, <<<nblk(Vx.entries()), TH, 0, s>>>(V.m, V.e, V.rows, V.cols, dim,
                                                  V.batch, Vx.m, Vx.e));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_repeat(const Mpf_Batch &A, Mpf_Batch &R, cudaStream_t s)
  {
    OPS_BY_L(A.L, k_repeat, <<<nblk(R.entries()), TH, 0, s>>>(
      A.m, A.e, A.entries(), R.m, R.e, int(R.batch / A.batch)));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_seg_diag(const Mpf_Batch &M, Mpf_Batch &v, int K, int dim,
                    cudaStream_t s)
  {
    OPS_BY_L(M.L, k_seg_diag, <<<nblk(v.entries()), TH, 0, s>>>(
      M.m, M.e, M.rows, M.batch / 2, K, dim * (dim + 1) / 2, v.m, v.e));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_scale_segments(const Mpf_Batch &VT, const Mpf_Batch &a, Mpf_Batch &VD,
                          int nseg, cudaStream_t s)
  {
    OPS_BY_L(VT.L, k_scale_segments, <<<nblk(VD.entries()), TH, 0, s>>>(
      VT.m, VT.e, VT.rows, VT.cols, VT.batch / 2, a.m, a.e, nseg, VD.m, VD.e));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_cmws_assemble(const Mpf_Batch &C, Mpf_Batch &R, int h, int dim,
                         cudaStream_t s)
  {
    OPS_BY_L(C.L, k_cmws, <<<nblk(R.entries()), TH, 0, s>>>(C.m, C.e, h, dim, R.batch, R.m,
                                                 R.e));
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_schur(const Mpf_Batch &AXI, const Mpf_Batch &AY, Mpf_Batch &S, int K,
                 int dim, const int *d_Kj, cudaStream_t s)
  {
    OPS_BY_L(AXI.L, k_schur, <<<nblk(S.entries()), TH, 0, s>>>(
      AXI.m, AXI.e, AY.m, AY.e, K, dim, S.batch, d_Kj, S.m, S.e));
    OPS_CHECK(cudaGetLastError());
  }

  // threads per output: about this many in all (SDPB_GPU_GEMV_THREADS; 0: one per output)
  int gemv_thread_target = [] {
    const char *e = std::getenv("SDPB_GPU_GEMV_THREADS");
    return e ? std::atoi(e) : 262144;
  }();

  void mpf_gemv(const Mpf_Batch &A, const Mpf_Batch &v, Mpf_Batch &out,
                const Gemv_Spec &g, cudaStream_t s)
  {
    const int64_t n = int64_t(g.batch) * g.no, nt = g.ntb * g.nti;
    int S = 1;
    if(gemv_thread_target > 0 && n > 0)
      S = int(std::max<int64_t>(1, std::min<int64_t>({64, gemv_thread_target / n,
                                                      nt / 256})));
    // few long-running threads spread over all the streaming multiprocessors
    int th = 32;
    while(th < TH && n / th > 132)
      th *= 2;
    const int64_t grid = (n + th - 1) / th;
    if(S > 1)
      {
        const int W = 2 * A.L + 2 - (A.L - 3);
        i64 *T;
        u32 *Part;
        OPS_CHECK(cudaMallocAsync(&T, size_t(n) * 8, s));
        OPS_CHECK(cudaMallocAsync(&Part, size_t(S) * W * n * 4, s));
        OPS_BY_L(A.L, k_gemv_top, <<<grid, th, 0, s>>>(A.e, v.e, g.no, g.batch, g.ab, g.ao,
                                            g.atb, g.ati, g.vb, g.vo, g.vtb, g.vti,
                                            g.ntb, g.nti, T));
        OPS_BY_L(A.L, k_gemv_part, <<<nblk(n * S), TH, 0, s>>>(
          A.m, A.e, A.entries(), v.m, v.e, v.entries(), T, g.no, g.batch, g.ab,
          g.ao, g.atb, g.ati, g.vb, g.vo, g.vtb, g.vti, g.ntb, g.nti, S,
          g.ati == 1 ? 1 : 0, Part));
        OPS_BY_L(A.L, k_gemv_reduce, <<<grid, th, 0, s>>>(Part, T, g.no, g.batch, S, g.ob,
                                               g.oo, out.m, out.e, out.entries()));
        OPS_CHECK(cudaGetLastError());
        OPS_CHECK(cudaFreeAsync(Part, s));
        OPS_CHECK(cudaFreeAsync(T, s));
        return;
      }
    OPS_BY_L(A.L, k_gemv, <<<grid, th, 0, s>>>(
      A.m, A.e, A.entries(), v.m, v.e, v.entries(), g.no, g.batch, g.ab, g.ao,
      g.atb, g.ati, g.vb, g.vo, g.vtb, g.vti, g.ntb, g.nti, g.ob, g.oo, out.m,
      out.e, out.entries()));
    OPS_CHECK(cudaGetLastError());
  }

  namespace
  {
    i64 max_exponent(const Mpf_Batch &A, cudaStream_t s)
    {
      u64 *d, h = 0;
      OPS_CHECK(cudaMallocAsync(&d, 8, s));
      OPS_CHECK(cudaMemsetAsync(d, 0, 8, s));
      k_max_exp2<<<TH, TH, 0, s>>>(A.e, A.entries(), d);
      OPS_CHECK(cudaGetLastError());
      OPS_CHECK(cudaMemcpyAsync(&h, d, 8, cudaMemcpyDeviceToHost, s));
      OPS_CHECK(cudaStreamSynchronize(s));
      OPS_CHECK(cudaFreeAsync(d, s));
      return (i64)(h ^ 0x8000000000000000ULL);
    }
  }

  Host_Float mpf_max_abs(const Mpf_Batch &A, cudaStream_t s)
  {
    check_L(A.L);
    Host_Float r;
    if(A.entries() == 0)
      return r;
    const i64 emax = max_exponent(A, s);
    if(emax == LLONG_MIN)
      return r;
    u64 *d, key = 0;
    OPS_CHECK(cudaMallocAsync(&d, 8, s));
    OPS_CHECK(cudaMemsetAsync(d, 0, 8, s));
    const int64_t n = A.entries();
    OPS_BY_L(A.L, k_max_key, <<<nblk(n), TH, 0, s>>>(A.m, A.e, n, emax, d));
    OPS_CHECK(cudaGetLastError());
    OPS_CHECK(cudaMemcpyAsync(&key, d, 8, cudaMemcpyDeviceToHost, s));
    OPS_CHECK(cudaStreamSynchronize(s));
    unsigned int cap = 256;
    int32_t *o;
    unsigned int *dcount, hcount = 0;
    OPS_CHECK(cudaMallocAsync(&o, 4 * size_t(A.L) * cap, s));
    OPS_CHECK(cudaMallocAsync(&dcount, 4, s));
    OPS_CHECK(cudaMemsetAsync(dcount, 0, 4, s));
    OPS_BY_L(A.L, k_find_key, <<<nblk(n), TH, 0, s>>>(A.m, A.e, n, emax, key, o, dcount, cap));
    OPS_CHECK(cudaGetLastError());
    OPS_CHECK(cudaMemcpyAsync(&hcount, dcount, 4, cudaMemcpyDeviceToHost, s));
    OPS_CHECK(cudaStreamSynchronize(s));
    // More ties in the top 64 bits than the buffer holds: which ones got a slot
    // depends on the threads, so the low digits of the maximum would depend on
    // the run. Now the count is known: all of them.
    if(hcount > cap)
      {
        cap = hcount;
        OPS_CHECK(cudaFreeAsync(o, s));
        OPS_CHECK(cudaMallocAsync(&o, 4 * size_t(A.L) * cap, s));
        OPS_CHECK(cudaMemsetAsync(dcount, 0, 4, s));
        OPS_BY_L(A.L, k_find_key, <<<nblk(n), TH, 0, s>>>(A.m, A.e, n, emax, key, o, dcount, cap));
        OPS_CHECK(cudaGetLastError());
        OPS_CHECK(cudaStreamSynchronize(s));
        if(const char *dbg = getenv("SDPB_GPU_DEBUG_TIES"); dbg != nullptr)
          fprintf(stderr, "mpf_max_abs: %u ties in the top 64 bits\n", hcount);
      }
    const unsigned int found = hcount;
    std::vector<uint32_t> buf(size_t(A.L) * found);
    OPS_CHECK(cudaMemcpyAsync(buf.data(), o, buf.size() * 4,
                              cudaMemcpyDeviceToHost, s));
    OPS_CHECK(cudaStreamSynchronize(s));
    size_t best = 0;
    for(size_t q = 1; q < found; ++q)
      for(int t = A.L - 1; t >= 0; --t)
        {
          const uint32_t x = buf[q * A.L + t], y = buf[best * A.L + t];
          if(x != y)
            {
              if(x > y)
                best = q;
              break;
            }
        }
    r.m.assign(buf.begin() + best * A.L, buf.begin() + (best + 1) * A.L);
    r.e = emax;
    OPS_CHECK(cudaFreeAsync(d, s));
    OPS_CHECK(cudaFreeAsync(o, s));
    OPS_CHECK(cudaFreeAsync(dcount, s));
    return r;
  }

  Host_Float mpf_total(const Mpf_Batch &A, cudaStream_t s)
  {
    check_L(A.L);
    Host_Float r;
    if(A.entries() == 0)
      return r;
    const i64 emax = max_exponent(A, s);
    if(emax == LLONG_MIN)
      return r;
    const i64 t0 = emax - 32;
    const int W = A.L + 2;
    u64 *d;
    std::vector<u64> h(W + 1, 0);
    OPS_CHECK(cudaMallocAsync(&d, 8 * (W + 1), s));
    OPS_CHECK(cudaMemsetAsync(d, 0, 8 * (W + 1), s));
    const int64_t n = A.entries();
    OPS_BY_L(A.L, k_total, <<<nblk(n), TH, 0, s>>>(A.m, A.e, n, t0, d));
    OPS_CHECK(cudaGetLastError());
    OPS_CHECK(cudaMemcpyAsync(h.data(), d, 8 * (W + 1), cudaMemcpyDeviceToHost, s));
    OPS_CHECK(cudaStreamSynchronize(s));
    OPS_CHECK(cudaFreeAsync(d, s));
    // sum_k h[k] 2^(32 k) as unsigned parts of two's complement W-limb terms;
    // sum over negative terms of -2^(32 W), h[W] counts them
    const int WW = W + 3;
    r.m.assign(WW, 0);
    u64 carry = 0;
    for(int k = 0; k < WW; ++k)
      {
        const u64 part = k < W ? h[k] : 0;
        const unsigned __int128 v = (unsigned __int128)part + carry;
        r.m[k] = u32(v);
        carry = u64(v >> 32);
      }
    // subtract h[W] 2^(32 W)
    unsigned __int128 borrow = (unsigned __int128)h[W];
    for(int k = W; k < WW; ++k)
      {
        const u32 sub = u32(borrow);
        const u32 old = r.m[k];
        r.m[k] = old - sub;
        borrow = (borrow >> 32) + (old < sub ? 1 : 0);
      }
    r.e = t0;
    return r;
  }

  void mpf_to_double(const Mpf_Batch &A, double *d_out, int64_t *scale,
                     cudaStream_t s)
  {
    check_L(A.L);
    const int64_t per = A.rows * A.cols;
    u64 *d;
    int64_t *dt;
    std::vector<u64> h(A.batch);
    OPS_CHECK(cudaMallocAsync(&d, 8 * size_t(A.batch), s));
    OPS_CHECK(cudaMallocAsync(&dt, 8 * size_t(A.batch), s));
    OPS_CHECK(cudaMemsetAsync(d, 0, 8 * size_t(A.batch), s));
    const unsigned gx = unsigned(std::max<int64_t>(1, (per + 4095) / 4096));
    k_block_top<<<dim3(gx, unsigned(A.batch)), TH, 0, s>>>(A.e, per, A.batch, d);
    OPS_CHECK(cudaGetLastError());
    OPS_CHECK(cudaMemcpyAsync(h.data(), d, 8 * size_t(A.batch),
                              cudaMemcpyDeviceToHost, s));
    OPS_CHECK(cudaStreamSynchronize(s));
    std::vector<int64_t> tops(A.batch);
    for(int b = 0; b < A.batch; ++b)
      {
        const i64 e = (i64)(h[b] ^ 0x8000000000000000ULL);
        // e + 32 (L - 2) + 62 bits of the two top limbs = e + NB
        tops[b] = (h[b] == 0 || e == LLONG_MIN) ? 0 : e + 32 * A.L - 2;
        scale[b] = tops[b];
      }
    OPS_CHECK(cudaMemcpyAsync(dt, tops.data(), 8 * size_t(A.batch),
                              cudaMemcpyHostToDevice, s));
    OPS_CHECK(cudaStreamSynchronize(s));
    const int64_t n = A.entries();
    OPS_BY_L(A.L, k_to_double, <<<nblk(n), TH, 0, s>>>(A.m, A.e, per, A.batch, dt, d_out));
    OPS_CHECK(cudaGetLastError());
    OPS_CHECK(cudaFreeAsync(d, s));
    OPS_CHECK(cudaFreeAsync(dt, s));
  }

  namespace
  {
    __global__ void k_unpack_sparse(const int32_t *vm, const int64_t *ve,
                                    const uint32_t *pos, int64_t nnz, int64_t tot,
                                    int L, int32_t *m, int64_t *e)
    {
      const int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(i >= nnz)
        return;
      const int64_t q = pos[i];
      e[q] = ve[i];
      for(int t = 0; t < L; ++t)
        m[int64_t(t) * tot + q] = vm[int64_t(t) * nnz + i];
    }
  }

  void mpf_pack_sparse(const int32_t *m, const int64_t *e, int64_t entries, int L,
                       std::vector<uint32_t> &pos, std::vector<int32_t> &vm,
                       std::vector<int64_t> &ve)
  {
    if(entries >= (int64_t(1) << 32))
      throw std::runtime_error("mpf_pack_sparse: more entries than 32-bit positions");
    pos.clear();
    for(int64_t q = 0; q < entries; ++q)
      {
        bool zero = e[q] == MPF_EZ;
        for(int t = 0; zero && t < L; ++t)
          zero = m[int64_t(t) * entries + q] == 0;
        if(!zero)
          pos.push_back(uint32_t(q));
      }
    const int64_t nnz = int64_t(pos.size());
    vm.resize(size_t(L) * nnz);
    ve.resize(nnz);
    for(int64_t i = 0; i < nnz; ++i)
      ve[i] = e[pos[i]];
    for(int t = 0; t < L; ++t)
      for(int64_t i = 0; i < nnz; ++i)
        vm[int64_t(t) * nnz + i] = m[int64_t(t) * entries + pos[i]];
  }

  void mpf_unpack_sparse(const Mpf_Batch &v, const uint32_t *pos, Mpf_Batch &A,
                         cudaStream_t s)
  {
    mpf_zero(A, s);
    const int64_t nnz = v.entries();
    const int th = 256;
    k_unpack_sparse<<<std::max<int64_t>(1, (nnz + th - 1) / th), th, 0, s>>>(
      v.m, v.e, pos, nnz, A.entries(), A.L, A.m, A.e);
    OPS_CHECK(cudaGetLastError());
  }

  void mpf_download(const Mpf_Batch &A, int32_t *m, int64_t *e, cudaStream_t s)
  {
    OPS_CHECK(cudaMemcpyAsync(m, A.m, size_t(A.L) * A.entries() * 4,
                              cudaMemcpyDeviceToHost, s));
    OPS_CHECK(cudaMemcpyAsync(e, A.e, size_t(A.entries()) * 8,
                              cudaMemcpyDeviceToHost, s));
    OPS_CHECK(cudaStreamSynchronize(s));
  }

  void mpf_upload(const int32_t *m, const int64_t *e, Mpf_Batch &A,
                  cudaStream_t s)
  {
    OPS_CHECK(cudaMemcpyAsync(A.m, m, size_t(A.L) * A.entries() * 4,
                              cudaMemcpyHostToDevice, s));
    OPS_CHECK(cudaMemcpyAsync(A.e, e, size_t(A.entries()) * 8,
                              cudaMemcpyHostToDevice, s));
    OPS_CHECK(cudaStreamSynchronize(s));
  }
}
