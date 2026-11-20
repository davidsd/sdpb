// GPU backend for SDPB: CUDA-core base case of the batched block
// factorizations: multiprecision Cholesky and triangular inverse of BT x BT
// tiles in wide fixed point, one thread block per tile.
//
// Ported from SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: sdpvibe/basecase.py (_SRC: mulfx, rsqrt_fx,
// the warp-cooperative wmul/wrsqrt, chol_tiles, trinv_tiles).  The NVRTC
// source with WL/LF/BT macros becomes C++ templates; the arithmetic is
// unchanged.  Portions adapted from SDPviBe, Copyright (c) 2026 Antoine Tilloy,
// MIT License (docs/GPU_backend_credits.md).
//
// Format: WL 32-bit limbs, two's complement, little endian, F = 32 LF
// fractional bits (LI = WL - LF integer limbs).  Tiles are entry-major:
// limb w of entry (i, j) of tile t at ((t*BT + i)*BT + j)*WL + w.
// Input tiles are SPD and equilibrated (|a_ij| <= 1); pivots may be tiny, and
// their reciprocals and the inverse X = L^-1 large (the integer limbs).
#pragma once

#include <cstdint>

namespace sdpb_gpu
{
  template <int WL, int LF, int BT> struct Tile_Kernels
  {
    typedef unsigned int u32;
    typedef unsigned long long u64;
    typedef long long i64;

    static __device__ __forceinline__ void cpy(u32 *r, const u32 *a)
    {
      for(int i = 0; i < WL; i++)
        r[i] = a[i];
    }
    static __device__ __forceinline__ void zero(u32 *r)
    {
      for(int i = 0; i < WL; i++)
        r[i] = 0;
    }
    static __device__ __forceinline__ bool isneg(const u32 *a)
    {
      return a[WL - 1] >> 31;
    }
    static __device__ void negate(u32 *a)
    {
      u64 c = 1;
      for(int i = 0; i < WL; i++)
        {
          u64 t = (u64)(~a[i]) + c;
          a[i] = (u32)t;
          c = t >> 32;
        }
    }
    static __device__ void addto(u32 *r, const u32 *a)
    {
      u64 c = 0;
      for(int i = 0; i < WL; i++)
        {
          u64 t = (u64)r[i] + a[i] + c;
          r[i] = (u32)t;
          c = t >> 32;
        }
    }
    static __device__ void subfrom(u32 *r, const u32 *a)
    {
      i64 br = 0;
      for(int i = 0; i < WL; i++)
        {
          i64 t = (i64)r[i] - (i64)a[i] + br;
          r[i] = (u32)t;
          br = t >> 32;
        }
    }
    static __device__ __forceinline__ int topl(const u32 *a)
    {
      for(int i = WL - 1; i > 0; i--)
        if(a[i])
          return i;
      return 0;
    }
    static __device__ void shr1(u32 *a)
    {
      for(int i = 0; i < WL - 1; i++)
        a[i] = (a[i] >> 1) | (a[i + 1] << 31);
      a[WL - 1] = (u32)(((int)a[WL - 1]) >> 1);
    }

    // r = (a*b) >> F, signed, truncated; t > 0: only the top t limbs of |a|, |b|
    static __device__ void mulfx(u32 *r, const u32 *a, const u32 *b, int t)
    {
      u32 x[WL], y[WL];
      cpy(x, a);
      cpy(y, b);
      bool sa = isneg(x), sb = isneg(y);
      if(sa)
        negate(x);
      if(sb)
        negate(y);
      int ta = topl(x), tb = topl(y);
      int la = t ? max(0, ta - t + 1) : 0, lb = t ? max(0, tb - t + 1) : 0;
      int kmin = max(LF - 2, la + lb), kmax = ta + tb;
      for(int j = 0; j < WL; j++)
        r[j] = 0;
      u64 carry = 0, hprev = 0;
      for(int k = kmin; k <= kmax + 2 && k < LF + WL; k++)
        {
          u64 lo = 0, hi = 0;
          if(k <= kmax)
            {
              int i0 = max(la, k - tb), i1 = min(ta, k - lb);
#pragma unroll 4
              for(int i = i0; i <= i1; i++)
                {
                  u64 p = (u64)x[i] * y[k - i];
                  lo += (u32)p;
                  hi += p >> 32;
                }
            }
          u64 s = lo + hprev + carry;
          if(k >= LF)
            r[k - LF] = (u32)s;
          carry = s >> 32;
          hprev = hi;
        }
      if(sa ^ sb)
        negate(r);
    }

    static __device__ double to_double(const u32 *a)
    {
      u32 x[WL];
      cpy(x, a);
      bool s = isneg(x);
      if(s)
        negate(x);
      int t = topl(x);
      double v = (double)x[t];
      if(t >= 1)
        v = v * 4294967296.0 + (double)x[t - 1];
      if(t >= 2)
        v = v * 4294967296.0 + (double)x[t - 2];
      int e = 32 * (max(t - 2, 0)) - 32 * LF;
      v = ldexp(v, e);
      return s ? -v : v;
    }
    static __device__ void from_double(u32 *r, double v) // v > 0
    {
      zero(r);
      int e;
      double m = frexp(v, &e);
      u64 mi = (u64)ldexp(m, 53);
      int pos = e - 53 + 32 * LF;
      if(pos < 0)
        {
          mi >>= (-pos);
          pos = 0;
        }
      int li = pos / 32, sh = pos % 32;
      unsigned __int128 w = ((unsigned __int128)mi) << sh;
      for(int k = 0; k < 4 && li + k < WL; k++)
        {
          r[li + k] = (u32)(w & 0xFFFFFFFFu);
          w >>= 32;
        }
    }

    // ---- warp-cooperative arithmetic (32 lanes; operands in shared memory)
    static __device__ __forceinline__ void wcpy(u32 *d, const u32 *a)
    {
      for(int i = threadIdx.x & 31; i < WL; i += 32)
        d[i] = a[i];
      __syncwarp();
    }
    static __device__ void wneg(u32 *a)
    {
      for(int i = threadIdx.x & 31; i < WL; i += 32)
        a[i] = ~a[i];
      __syncwarp();
      if((threadIdx.x & 31) == 0)
        {
          for(int i = 0; i < WL; i++)
            {
              a[i] += 1u;
              if(a[i])
                break;
            }
        }
      __syncwarp();
    }
    static __device__ void wmul(u32 *r, const u32 *a, const u32 *b, int t,
                                u32 *sx, u32 *sy, u64 *col)
    {
      int lane = threadIdx.x & 31;
      wcpy(sx, a);
      wcpy(sy, b);
      bool sa = isneg(sx), sb = isneg(sy);
      __syncwarp();
      if(sa)
        wneg(sx);
      if(sb)
        wneg(sy);
      int ta = topl(sx), tb = topl(sy);
      int la = t ? max(0, ta - t + 1) : 0, lb = t ? max(0, tb - t + 1) : 0;
      int kmin = max(LF - 2, la + lb), kmax = ta + tb;
      for(int k = kmin + lane; k <= kmax; k += 32)
        {
          u64 lo = 0, hi = 0;
          int i0 = max(la, k - tb), i1 = min(ta, k - lb);
          for(int i = i0; i <= i1; i++)
            {
              u64 p = (u64)sx[i] * sy[k - i];
              lo += (u32)p;
              hi += p >> 32;
            }
          col[2 * (k - kmin)] = lo;
          col[2 * (k - kmin) + 1] = hi;
        }
      for(int j = lane; j < WL; j += 32)
        r[j] = 0;
      __syncwarp();
      if(lane == 0)
        {
          u64 carry = 0, hprev = 0;
          for(int k = kmin; k <= kmax + 2 && k < LF + WL; k++)
            {
              u64 lo = 0, hi = 0;
              if(k <= kmax)
                {
                  lo = col[2 * (k - kmin)];
                  hi = col[2 * (k - kmin) + 1];
                }
              u64 s = lo + hprev + carry;
              if(k >= LF)
                r[k - LF] = (u32)s;
              carry = s >> 32;
              hprev = hi;
            }
        }
      __syncwarp();
      if(sa ^ sb)
        wneg(r);
    }
    static __device__ void wsub(u32 *r, const u32 *a)
    {
      if((threadIdx.x & 31) == 0)
        subfrom(r, a);
      __syncwarp();
    }
    static __device__ void wadd(u32 *r, const u32 *a)
    {
      if((threadIdx.x & 31) == 0)
        addto(r, a);
      __syncwarp();
    }
    static __device__ void wshr1(u32 *a)
    {
      int lane = threadIdx.x & 31;
      u32 v[(WL + 31) / 32];
      int c = 0;
      for(int i = lane; i < WL; i += 32, c++)
        v[c] = (i < WL - 1) ? ((a[i] >> 1) | (a[i + 1] << 31))
                            : (u32)(((int)a[i]) >> 1);
      __syncwarp();
      c = 0;
      for(int i = lane; i < WL; i += 32, c++)
        a[i] = v[c];
      __syncwarp();
    }
    // r = d^{-1/2}, precision-doubling Newton r <- r + r (1 - d r^2) / 2
    static __device__ void wrsqrt(u32 *r, const u32 *d, u32 *t1, u32 *t2,
                                  u32 *one, u32 *sx, u32 *sy, u64 *col)
    {
      int lane = threadIdx.x & 31;
      if(lane == 0)
        {
          from_double(r, 1.0 / sqrt(to_double(d)));
          zero(one);
          one[LF] = 1;
        }
      __syncwarp();
      for(int bits = 100;; bits *= 2)
        {
          int t = bits / 32 + 2;
          if(t >= WL)
            t = 0;
          wmul(t1, r, r, t, sx, sy, col);
          wmul(t2, t1, d, t, sx, sy, col);
          wcpy(t1, one);
          wsub(t1, t2);
          wmul(t2, r, t1, t, sx, sy, col);
          wshr1(t2);
          wadd(r, t2);
          if(t == 0)
            break;
        }
    }
  };

  // batched Cholesky of BT x BT tiles (in place, lower) + pivot reciprocals
  template <int WL, int LF, int BT>
  __global__ void chol_tiles(unsigned int *A, unsigned int *RINV)
  {
    using K = Tile_Kernels<WL, LF, BT>;
    typedef unsigned int u32;
    typedef unsigned long long u64;
    int tile = blockIdx.x;
    u32 *T = A + (size_t)tile * BT * BT * WL;
    u32 *RI = RINV + (size_t)tile * BT * WL;
    __shared__ u32 sd[WL], sr[WL], st1[WL], st2[WL], sone[WL], sx[WL], sy[WL];
    __shared__ u64 scol[2 * (2 * WL + 4)];
    for(int k = 0; k < BT; k++)
      {
        if(threadIdx.x < 32)
          {
            u32 *akk = T + (size_t)(k * BT + k) * WL;
            K::wcpy(sd, akk);
            K::wrsqrt(sr, sd, st1, st2, sone, sx, sy, scol);
            K::wmul(st1, sd, sr, 0, sx, sy, scol);
            K::wcpy(akk, st1);
            K::wcpy(RI + (size_t)k * WL, sr);
          }
        __syncthreads();
        const u32 *rk = RI + (size_t)k * WL;
        for(int i = k + 1 + threadIdx.x; i < BT; i += blockDim.x)
          {
            u32 *aik = T + (size_t)(i * BT + k) * WL;
            u32 v[WL];
            K::mulfx(v, aik, rk, 0);
            K::cpy(aik, v);
          }
        __syncthreads();
        int m = BT - k - 1, npairs = m * (m + 1) / 2;
        for(int p = threadIdx.x; p < npairs; p += blockDim.x)
          {
            int ii = (int)((sqrtf(8.0f * p + 1.0f) - 1.0f) * 0.5f);
            while((ii + 1) * (ii + 2) / 2 <= p)
              ii++;
            while(ii * (ii + 1) / 2 > p)
              ii--;
            int jj = p - ii * (ii + 1) / 2;
            int i = k + 1 + ii, j = k + 1 + jj;
            u32 v[WL];
            K::mulfx(v, T + (size_t)(i * BT + k) * WL,
                     T + (size_t)(j * BT + k) * WL, 0);
            K::subfrom(T + (size_t)(i * BT + j) * WL, v);
          }
        __syncthreads();
      }
  }

  // batched X = L^-1 (lower) from the factor and the pivot reciprocals:
  // x_ii = 1/l_ii, x_ik = -(1/l_ii) sum_{j=k}^{i-1} l_ij x_jk
  template <int WL, int LF, int BT>
  __global__ void trinv_tiles(const unsigned int *Lt, const unsigned int *RINV,
                              unsigned int *X, unsigned int *SCR)
  {
    using K = Tile_Kernels<WL, LF, BT>;
    typedef unsigned int u32;
    int tile = blockIdx.x;
    const u32 *T = Lt + (size_t)tile * BT * BT * WL;
    const u32 *RI = RINV + (size_t)tile * BT * WL;
    u32 *XT = X + (size_t)tile * BT * BT * WL;
    u32 *S = SCR + (size_t)tile * BT * BT * WL;
    for(int e = threadIdx.x; e < BT * BT; e += blockDim.x)
      K::zero(XT + (size_t)e * WL);
    __syncthreads();
    for(int i = 0; i < BT; i++)
      {
        int np = i * (i + 1) / 2;
        for(int p = threadIdx.x; p < np; p += blockDim.x)
          {
            int jj = (int)((sqrtf(8.0f * p + 1.0f) - 1.0f) * 0.5f);
            while((jj + 1) * (jj + 2) / 2 <= p)
              jj++;
            while(jj * (jj + 1) / 2 > p)
              jj--;
            int k = p - jj * (jj + 1) / 2, j = jj;
            K::mulfx(S + (size_t)(k * BT + j) * WL, T + (size_t)(i * BT + j) * WL,
                     XT + (size_t)(j * BT + k) * WL, 0);
          }
        __syncthreads();
        for(int k = threadIdx.x; k <= i; k += blockDim.x)
          {
            u32 *xik = XT + (size_t)(i * BT + k) * WL;
            if(k == i)
              {
                K::cpy(xik, RI + (size_t)i * WL);
                continue;
              }
            u32 acc[WL];
            K::zero(acc);
            for(int j = k; j < i; j++)
              K::addto(acc, S + (size_t)(k * BT + j) * WL);
            u32 v[WL];
            K::mulfx(v, acc, RI + (size_t)i * WL, 0);
            K::negate(v);
            K::cpy(xik, v);
          }
        __syncthreads();
      }
  }
}
