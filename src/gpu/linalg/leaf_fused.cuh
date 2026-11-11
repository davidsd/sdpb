// GPU backend for SDPB: fused base case of the batched block factorization:
// the Cholesky factor L of a BT x BT tile and its inverse X = L^-1 in one
// right-looking pass, one thread block per tile.
//
// Same number format and conventions as basecase.cuh (SDPviBe's base case,
// ported from SDPviBe by Antoine Tilloy, https://github.com/AntoineTilloy/SDPviBe,
// sdpvibe/basecase.py: WL limbs of wide fixed point with F = 32 LF fractional
// bits, products truncated below column LF - 2, the pivot reciprocal square
// root by precision-doubling Newton from a double); the pivot code below uses
// the ported Tile_Kernels as is.  Portions adapted from SDPviBe, Copyright (c)
// 2026 Antoine Tilloy, MIT License (docs/GPU_backend_credits.md).
//
// New here: SDPviBe's base case factors a tile in BT
// sequential column steps and then inverts the factor in BT more row steps,
// with the products in each step spread over 256 threads that copy their
// operands to local memory.  Here the inverse is accumulated column by
// column alongside the factorization (forward substitution L X = I, rows of X
// finished as the pivots appear), so a tile takes BT steps instead of 2 BT;
// each step's products are one per thread (up to 527 for BT = 32; with 512
// threads, steps 0-4 take a second round for the last <= 15), their
// operands read from shared memory (sign-magnitude, with top-limb bounds),
// their column sums in a three-word carry chain and their output limbs in
// registers.
//
// Step k (r = a_kk^-1/2, so that 1/l_kk = r):
//   A (one warp)    r by Newton; l_kk = a_kk r
//   B (BT groups)   Y_s = r x_{k,s} accumulated (s < k), Y_k = r = x_kk,
//                   Y_s = r a_sk = l_sk (s > k); written to X and L
//   C (all)         for i > k, 0 <= s <= i: T_is -= Y_i Y_s, with T = X
//                   (the accumulated row i of the inverse) for s <= k and T = A
//                   (the trailing Schur complement) for s > k
#pragma once

#include "basecase.cuh"

namespace sdpb_gpu
{
  namespace leaf_detail
  {
    typedef unsigned int u32;

    // (c2:c1:c0) += a b
    __device__ __forceinline__ void mac3(u32 &c0, u32 &c1, u32 &c2, u32 a,
                                         u32 b)
    {
      asm("mad.lo.cc.u32 %0, %3, %4, %0;\n\t"
          "madc.hi.cc.u32 %1, %3, %4, %1;\n\t"
          "addc.u32 %2, %2, 0;"
          : "+r"(c0), "+r"(c1), "+r"(c2)
          : "r"(a), "r"(b));
    }

    // res = (x y) >> F for magnitudes x, y (top nonzero limbs ta, tb; -1 if
    // zero), columns below LF - 2 dropped as in Tile_Kernels::mulfx (t = 0)
    template <int WL, int LF>
    __device__ __forceinline__ void mag_mul(u32 (&res)[WL], const u32 *x,
                                            int ta, const u32 *y, int tb)
    {
      u32 c0 = 0, c1 = 0, c2 = 0;
#pragma unroll
      for(int k = LF - 2; k < LF + WL; ++k)
        {
          const int i0 = max(0, k - tb), i1 = min(ta, k);
          // two independent chains for latency, merged at the column's end
          u32 d0 = 0, d1 = 0, d2 = 0;
          int i = i0;
          for(; i + 1 <= i1; i += 2)
            {
              const u32 xa = x[i], ya = y[k - i], xb = x[i + 1],
                        yb = y[k - i - 1];
              mac3(c0, c1, c2, xa, ya);
              mac3(d0, d1, d2, xb, yb);
            }
          if(i <= i1)
            mac3(c0, c1, c2, x[i], y[k - i]);
          asm("add.cc.u32 %0, %0, %3;\n\t"
              "addc.cc.u32 %1, %1, %4;\n\t"
              "addc.u32 %2, %2, %5;"
              : "+r"(c0), "+r"(c1), "+r"(c2)
              : "r"(d0), "r"(d1), "r"(d2));
          if(k >= LF)
            res[k - LF] = c0;
          c0 = c1;
          c1 = c2;
          c2 = 0;
        }
    }

    // the same for ta < NX, tb < NY (x an entry of L, |x| < 2^32): x in
    // registers, a fixed set of products (zero limbs contribute nothing), no
    // data-dependent branches
    template <int WL, int LF, int NX, int NY = NX>
    __device__ __forceinline__ void mag_mul_small(u32 (&res)[WL], const u32 *x,
                                                  const u32 *y)
    {
      static_assert(NX - 1 + NY - 1 < LF + WL, "NX + NY");
      u32 xr[NX];
#pragma unroll
      for(int i = 0; i < NX; ++i)
        xr[i] = x[i];
      u32 c0 = 0, c1 = 0, c2 = 0;
#pragma unroll
      for(int k = LF - 2; k < LF + WL; ++k)
        {
#pragma unroll
          for(int i = (k - (NY - 1) > 0 ? k - (NY - 1) : 0);
              i <= (k < NX - 1 ? k : NX - 1); ++i)
            mac3(c0, c1, c2, xr[i], y[k - i]);
          if(k >= LF)
            res[k - LF] = c0;
          c0 = c1;
          c1 = c2;
          c2 = 0;
        }
    }

    template <int WL>
    __device__ __forceinline__ int top_reg(const u32 (&a)[WL])
    {
      int t = -1;
#pragma unroll
      for(int w = 0; w < WL; ++w)
        t = a[w] ? w : t;
      return t;
    }

    // t (global, two's complement) -= v (add: += v), v a magnitude; four
    // limbs at a time
    template <int WL>
    __device__ __forceinline__ void rmw(u32 *t, const u32 (&v)[WL], bool add)
    {
      static_assert(WL % 4 == 0, "WL % 4");
      uint4 *t4 = reinterpret_cast<uint4 *>(t);
      long long c = 0;
      auto step = [&](u32 &a, u32 b) {
        const long long s = add ? (long long)a + b + c : (long long)a - b + c;
        a = u32(s);
        c = s >> 32;
      };
#pragma unroll
      for(int q = 0; q < WL / 4; ++q)
        {
          uint4 z = t4[q];
          step(z.x, v[4 * q]);
          step(z.y, v[4 * q + 1]);
          step(z.z, v[4 * q + 2]);
          step(z.w, v[4 * q + 3]);
          t4[q] = z;
        }
    }

    // global two's complement -> shared magnitude, sign, top limb
    template <int WL>
    __device__ __forceinline__ void load_mag(const u32 *g, u32 *m, int &top,
                                             int &neg)
    {
      u32 a[WL];
      const uint4 *g4 = reinterpret_cast<const uint4 *>(g);
#pragma unroll
      for(int q = 0; q < WL / 4; ++q)
        {
          const uint4 z = g4[q];
          a[4 * q] = z.x;
          a[4 * q + 1] = z.y;
          a[4 * q + 2] = z.z;
          a[4 * q + 3] = z.w;
        }
      neg = int(a[WL - 1] >> 31);
      unsigned long long c = 1;
#pragma unroll
      for(int w = 0; w < WL; ++w)
        if(neg)
          {
            const unsigned long long s = (unsigned long long)(~a[w]) + c;
            a[w] = u32(s);
            c = s >> 32;
          }
#pragma unroll
      for(int w = 0; w < WL; ++w)
        m[w] = a[w];
      top = top_reg<WL>(a);
    }

    // magnitude (registers) with sign -> global two's complement
    template <int WL>
    __device__ __forceinline__ void store_signed(u32 *g, const u32 (&v)[WL],
                                                 bool neg)
    {
      uint4 *g4 = reinterpret_cast<uint4 *>(g);
      unsigned long long c = 1;
      auto f = [&](u32 b) -> u32 {
        if(!neg)
          return b;
        const unsigned long long s = (unsigned long long)(~b) + c;
        c = s >> 32;
        return u32(s);
      };
#pragma unroll
      for(int q = 0; q < WL / 4; ++q)
        {
          const u32 a0 = f(v[4 * q]), a1 = f(v[4 * q + 1]),
                    a2 = f(v[4 * q + 2]), a3 = f(v[4 * q + 3]);
          g4[q] = make_uint4(a0, a1, a2, a3);
        }
    }

    // a_kk -> r = a_kk^-1/2 (R: magnitude, rtop: its top limb) and, in
    // place, l_kk = a_kk r: Newton as in Tile_Kernels::wrsqrt, one thread
    template <int WL, int LF, int BT>
    __device__ __noinline__ void pivot_rsqrt(u32 *akk, u32 *R, int *rtop)
    {
      using K = Tile_Kernels<WL, LF, BT>;
      u32 d[WL], r[WL], t1[WL], t2[WL];
      K::cpy(d, akk);
      K::from_double(r, 1.0 / sqrt(K::to_double(d)));
      for(int bits = 100;; bits *= 2)
        {
          int t = bits / 32 + 2;
          if(t >= WL)
            t = 0;
          K::mulfx(t1, d, r, t);  // d r
          K::mulfx(t2, t1, r, t); // d r^2
          K::zero(t1);
          t1[LF] = 1;
          K::subfrom(t1, t2);     // 1 - d r^2
          K::mulfx(t2, r, t1, t);
          K::shr1(t2);
          K::addto(r, t2);        // r += r (1 - d r^2) / 2
          if(t == 0)
            break;
        }
      K::mulfx(t1, d, r, 0);      // l_kk = a_kk r
      K::cpy(akk, t1);
      K::cpy(R, r);
      *rtop = K::topl(r);
    }

    // ---- products and sums spread over a group of G lanes (G | 32, the group
    // G-aligned in its warp; every lane of the group calls with the same
    // arguments).  Operands and results are magnitudes in shared memory.
    template <int G> __device__ __forceinline__ unsigned group_mask()
    {
      return G == 32 ? 0xffffffffu
                     : (((1u << G) - 1) << ((threadIdx.x & 31) & ~(G - 1)));
    }
    template <int G> __device__ __forceinline__ int group_lane()
    {
      return threadIdx.x & (G - 1);
    }
    // value at position j - d (d = 1, 2) of a strided vector: position
    // j = g + G m lives in lane g, slot m
    template <int G, int M>
    __device__ __forceinline__ void shift_up(const u32 (&v)[M], u32 (&o)[M],
                                             int d, unsigned mask)
    {
      const int g = group_lane<G>(), src = (g - d) & (G - 1);
#pragma unroll
      for(int m = 0; m < M; ++m)
        {
          const u32 a = __shfl_sync(mask, v[m], src, G);
          const u32 b = m ? __shfl_sync(mask, v[m ? m - 1 : 0], src, G) : 0u;
          o[m] = g >= d ? a : b;
        }
    }
    // limbs (u32) plus small carries c (value at position j + 1): propagate
    template <int G, int M>
    __device__ __forceinline__ void carry_out(u32 (&limb)[M], u32 (&c)[M],
                                              unsigned mask)
    {
      for(;;)
        {
          bool any = false;
#pragma unroll
          for(int m = 0; m < M; ++m)
            any = any || c[m];
          if(!__any_sync(mask, any))
            break;
          u32 ci[M];
          shift_up<G, M>(c, ci, 1, mask);
#pragma unroll
          for(int m = 0; m < M; ++m)
            {
              const unsigned long long s = (unsigned long long)limb[m] + ci[m];
              limb[m] = u32(s);
              c[m] = u32(s >> 32);
            }
        }
    }
    // top nonzero limb of a strided vector (output limb p = j + off; limbs
    // at p < 0 are not part of the output), -1 if 0
    template <int G, int M>
    __device__ __forceinline__ int group_top(const u32 (&limb)[M], int off,
                                             unsigned mask)
    {
      unsigned t = 0;
#pragma unroll
      for(int m = 0; m < M; ++m)
        {
          const int p = group_lane<G>() + G * m + off;
          if(limb[m] && p >= 0)
            t = max(t, unsigned(p + 1));
        }
      return int(__reduce_max_sync(mask, t)) - 1;
    }

    // out = (x y) >> F, as mulfx(out, x, y, t): only the top t limbs of x and
    // y (all if t = 0), columns below max(LF - 2, la + lb) dropped.  Column k
    // in lane (k - kmin) mod G as a three-word sum, then the words moved up to
    // their limbs and the carries propagated across the group.  Returns the
    // top limb of out.
    template <int WL, int LF, int G>
    __device__ int gmag_mul(u32 *out, const u32 *x, int ta, const u32 *y,
                            int tb, int t)
    {
      constexpr int M = (WL + 4 + G - 1) / G;
      const unsigned mask = group_mask<G>();
      const int g = group_lane<G>();
      for(int w = g; w < WL; w += G)
        out[w] = 0;
      if(ta < 0 || tb < 0)
        {
          __syncwarp(mask);
          return -1;
        }
      const int la = t ? max(0, ta - t + 1) : 0, lb = t ? max(0, tb - t + 1) : 0;
      const int kmin = max(LF - 2, la + lb), kmax = min(ta + tb, LF + WL - 1);
      u32 w0[M], w1[M], w2[M];
#pragma unroll
      for(int m = 0; m < M; ++m)
        {
          const int k = kmin + g + G * m;
          u32 c0 = 0, c1 = 0, c2 = 0;
          if(k <= kmax)
            {
              const int i0 = max(la, k - tb), i1 = min(ta, k - lb);
              u32 d0 = 0, d1 = 0, d2 = 0;
              int i = i0;
              for(; i + 1 <= i1; i += 2)
                {
                  const u32 xa = x[i], ya = y[k - i], xb = x[i + 1],
                            yb = y[k - i - 1];
                  mac3(c0, c1, c2, xa, ya);
                  mac3(d0, d1, d2, xb, yb);
                }
              if(i <= i1)
                mac3(c0, c1, c2, x[i], y[k - i]);
              asm("add.cc.u32 %0, %0, %3;\n\t"
                  "addc.cc.u32 %1, %1, %4;\n\t"
                  "addc.u32 %2, %2, %5;"
                  : "+r"(c0), "+r"(c1), "+r"(c2)
                  : "r"(d0), "r"(d1), "r"(d2));
            }
          w0[m] = c0;
          w1[m] = c1;
          w2[m] = c2;
        }
      u32 p1[M], p2[M], limb[M], c[M];
      shift_up<G, M>(w1, p1, 1, mask);
      shift_up<G, M>(w2, p2, 2, mask);
#pragma unroll
      for(int m = 0; m < M; ++m)
        {
          const unsigned long long u
            = (unsigned long long)w0[m] + p1[m] + p2[m];
          limb[m] = u32(u);
          c[m] = u32(u >> 32);
        }
      carry_out<G, M>(limb, c, mask);
      __syncwarp(mask);
#pragma unroll
      for(int m = 0; m < M; ++m)
        {
          const int p = kmin + g + G * m - LF;
          if(p >= 0 && p < WL)
            out[p] = limb[m];
        }
      const int top = group_top<G, M>(limb, kmin - LF, mask);
      __syncwarp(mask);
      return min(top, WL - 1);
    }

    // one warp: R = a_kk^-1/2 by r <- r (3 - a_kk r^2) / 2 (the Newton step
    // of Tile_Kernels::wrsqrt, same precision schedule), l_kk = a_kk r written
    // to akk.  Dm, T1, T2: shared scratch of WL limbs.
    // loaded: a_kk is already in Dm (look-ahead), not read from akk
    template <int WL, int LF, int BT>
    __device__ void pivot_warp(u32 *akk, u32 *R, int *rtop, u32 *Dm, u32 *T1,
                               u32 *T2, bool loaded = false)
    {
      using K = Tile_Kernels<WL, LF, BT>;
      constexpr int M = (WL + 31) / 32;
      const unsigned mask = 0xffffffffu;
      const int g = threadIdx.x & 31;
      u32 v[M];
#pragma unroll
      for(int m = 0; m < M; ++m)
        {
          const int p = g + 32 * m;
          v[m] = p < WL ? (loaded ? Dm[p] : akk[p]) : 0u;
          if(p < WL)
            {
              Dm[p] = v[m];
              R[p] = 0;
            }
        }
      const int dtop = group_top<32, M>(v, 0, mask);
      __syncwarp();
      if(g == 0)
        {
          u32 r0[WL];
          K::from_double(r0, 1.0 / sqrt(K::to_double(Dm)));
          for(int p = 0; p < WL; ++p)
            if(r0[p])
              R[p] = r0[p];
        }
      __syncwarp();
      int top = 0;
      for(int p = WL - 1; p >= 0; --p)
        if(R[p])
          {
            top = p;
            break;
          }
      for(int bits = 100;; bits *= 2)
        {
          int t = bits / 32 + 2;
          if(t >= WL)
            t = 0;
          int t1 = gmag_mul<WL, LF, 32>(T1, Dm, dtop, R, top, t); // d r
          int t2 = gmag_mul<WL, LF, 32>(T2, T1, t1, R, top, t);   // d r^2
          // T1 = 3 - T2 (positive)
          u32 l[M], b[M];
#pragma unroll
          for(int m = 0; m < M; ++m)
            {
              const int p = g + 32 * m;
              const long long s
                = (long long)(p == LF ? 3u : 0u) - (p < WL ? T2[p] : 0u);
              l[m] = u32(s);
              b[m] = s < 0;
            }
          for(;;)
            {
              bool any = false;
#pragma unroll
              for(int m = 0; m < M; ++m)
                any = any || b[m];
              if(!__any_sync(mask, any))
                break;
              u32 bi[M];
              shift_up<32, M>(b, bi, 1, mask);
#pragma unroll
              for(int m = 0; m < M; ++m)
                {
                  const long long s = (long long)l[m] - bi[m];
                  l[m] = u32(s);
                  b[m] = s < 0;
                }
            }
#pragma unroll
          for(int m = 0; m < M; ++m)
            if(g + 32 * m < WL)
              T1[g + 32 * m] = l[m];
          t1 = group_top<32, M>(l, 0, mask);
          __syncwarp();
          gmag_mul<WL, LF, 32>(T2, R, top, T1, t1, t); // r (3 - d r^2)
          // R = T2 / 2
#pragma unroll
          for(int m = 0; m < M; ++m)
            {
              const int p = g + 32 * m;
              l[m] = p < WL ? T2[p] : 0u;
            }
#pragma unroll
          for(int m = 0; m < M; ++m)
            {
              const int p = g + 32 * m;
              const u32 hi_same = __shfl_down_sync(mask, l[m], 1);
              const u32 hi_next = __shfl_sync(mask, l[m + 1 < M ? m + 1 : m], 0);
              u32 hi = g < 31 ? hi_same : (m + 1 < M ? hi_next : 0u);
              if(p + 1 >= WL)
                hi = 0;
              if(p < WL)
                R[p] = (l[m] >> 1) | (hi << 31);
            }
          __syncwarp();
#pragma unroll
          for(int m = 0; m < M; ++m)
            {
              const int p = g + 32 * m;
              l[m] = p < WL ? R[p] : 0u;
            }
          top = group_top<32, M>(l, 0, mask);
          if(t == 0)
            break;
        }
      gmag_mul<WL, LF, 32>(T1, Dm, dtop, R, top, 0); // l_kk = d r
      for(int p = g; p < WL; p += 32)
        akk[p] = T1[p];
      if(g == 0)
        *rtop = top;
      __syncwarp();
    }

    // look-ahead (one warp): d = a_jj - Y_j^2 (the last update of the next
    // pivot, j = k + 1) into Dm, then pivot_warp on it
    template <int WL, int LF, int BT>
    __device__ void next_pivot(u32 *ajj, const u32 *Yj, int ytop_j, u32 *R,
                               int *rtop, u32 *Dm, u32 *T1, u32 *T2)
    {
      constexpr int M = (WL + 31) / 32;
      const unsigned mask = 0xffffffffu;
      const int g = threadIdx.x & 31;
      gmag_mul<WL, LF, 32>(T1, Yj, ytop_j, Yj, ytop_j, 0);
      u32 l[M], b[M];
#pragma unroll
      for(int m = 0; m < M; ++m)
        {
          const int p = g + 32 * m;
          const long long sv
            = p < WL ? (long long)ajj[p] - (long long)T1[p] : 0LL;
          l[m] = u32(sv);
          b[m] = sv < 0;
        }
      for(;;)
        {
          bool any = false;
#pragma unroll
          for(int m = 0; m < M; ++m)
            any = any || b[m];
          if(!__any_sync(mask, any))
            break;
          u32 bi[M];
          shift_up<32, M>(b, bi, 1, mask);
#pragma unroll
          for(int m = 0; m < M; ++m)
            {
              const long long sv = (long long)l[m] - bi[m];
              l[m] = u32(sv);
              b[m] = sv < 0;
            }
        }
#pragma unroll
      for(int m = 0; m < M; ++m)
        if(g + 32 * m < WL)
          Dm[g + 32 * m] = l[m];
      __syncwarp();
      pivot_warp<WL, LF, BT>(ajj, R, rtop, Dm, T1, T2, true);
    }

    // a group of G lanes: global two's complement -> shared magnitude; returns
    // the top limb, sets neg
    template <int WL, int G>
    __device__ int group_load_mag(const u32 *src, u32 *dst, bool &neg)
    {
      const unsigned mask = group_mask<G>();
      const int g = group_lane<G>();
      neg = int(src[WL - 1]) < 0;
      unsigned low = WL;
      for(int w = g; w < WL; w += G)
        if(src[w])
          low = min(low, unsigned(w));
      low = __reduce_min_sync(mask, low);
      unsigned top = 0;
      for(int w = g; w < WL; w += G)
        {
          const u32 v = src[w];
          const u32 mg = !neg ? v : unsigned(w) < low ? 0u : unsigned(w) == low ? 0u - v : ~v;
          dst[w] = mg;
          if(mg)
            top = max(top, unsigned(w + 1));
        }
      __syncwarp(mask);
      return int(__reduce_max_sync(mask, top)) - 1;
    }
  }

  // A: tiles (entry-major, as in basecase.cuh), SPD and equilibrated; on exit
  // its lower triangle holds L.  X: on exit the lower triangle of L^-1, upper
  // triangle zero.  prof (optional): per tile, clock cycles of thread 0 in
  // phases A, B, C summed over the steps.
  // WARP: the pivot by one warp and row/column k by groups of NT / BT lanes
  // (gmag_mul) instead of one thread each
  // LA (with WARP): look-ahead: warp 0 updates and factors the next pivot
  // during phase C, which the other warps do; phase B loads its own operands
  template <int WL, int LF, int BT, int NT, bool WARP = true, bool LA = WARP>
  __global__ void __launch_bounds__(NT, 1)
    leaf_chol_inv(unsigned int *A, unsigned int *X, long long *prof)
  {
    using namespace leaf_detail;
    using K = Tile_Kernels<WL, LF, BT>;
    static_assert(NT >= 2 * BT, "phase B loaders");
    constexpr int WS = WL + 1; // odd stride: rows of Y in distinct banks
    __shared__ u32 Y[BT][WS], V[BT][WS], R[WS];
    __shared__ u32 Dm[WS], T1[WS], T2[WS];
    __shared__ int ytop[BT], yneg[BT], vtop[BT], vneg[BT], rtop;
    constexpr int GB = NT / BT; // lanes per product in phase B (WARP)
    static_assert(!WARP || (GB >= 1 && GB <= 32 && (GB & (GB - 1)) == 0),
                  "NT / BT a power of two <= 32");

    const int tid = threadIdx.x;
    u32 *T = A + (size_t)blockIdx.x * BT * BT * WL;
    u32 *XT = X + (size_t)blockIdx.x * BT * BT * WL;
    auto at = [&](u32 *M, int i, int j) { return M + (size_t)(i * BT + j) * WL; };

    for(int e = tid; e < BT * BT; e += NT)
      {
        u32 *x = XT + (size_t)e * WL;
        for(int w = 0; w < WL; ++w)
          x[w] = 0;
        if(e / BT == e % BT)
          x[LF] = 1;
      }
    __syncthreads();

    static_assert(!LA || WARP, "look-ahead needs WARP");
    long long cyc[3] = {0, 0, 0}, c0 = 0;
    if(LA)
      {
        if(prof && tid == 0)
          c0 = clock64();
        if(tid < 32)
          pivot_warp<WL, LF, BT>(at(T, 0, 0), R, &rtop, Dm, T1, T2);
        __syncthreads();
        if(prof && tid == 0)
          cyc[0] += clock64() - c0;
      }
    for(int k = 0; k < BT; ++k)
      {
        if(prof && tid == 0)
          c0 = clock64();
        // ---- A: r = a_kk^-1/2 (thread 0); operands of B loaded meanwhile
        if(LA)
          ;
        else if(WARP && tid < 32)
          pivot_warp<WL, LF, BT>(at(T, k, k), R, &rtop, Dm, T1, T2);
        else if(!WARP && tid == 0)
          pivot_rsqrt<WL, LF, BT>(at(T, k, k), R, &rtop);
        else if(tid >= 32 && tid < 32 + BT)
          {
            const int s = tid - 32;
            if(s != k)
              load_mag<WL>(s < k ? at(XT, k, s) : at(T, s, k), V[s], vtop[s],
                           vneg[s]);
          }
        if(!LA)
          __syncthreads();
        if(prof && tid == 0)
          {
            const long long c1 = clock64();
            cyc[0] += c1 - c0;
            c0 = c1;
          }
        // ---- B: row k of X, column k of L
        if(WARP)
          {
            const int s = tid / GB, g = tid % GB;
            if(LA && s != k)
              {
                bool ng;
                const int top = group_load_mag<WL, GB>(
                  s < k ? at(XT, k, s) : at(T, s, k), V[s], ng);
                if(g == 0)
                  {
                    vtop[s] = top;
                    vneg[s] = ng;
                  }
                __syncwarp(group_mask<GB>());
              }
            if(s == k)
              {
                for(int w = g; w < WL; w += GB)
                  Y[s][w] = R[w];
                if(g == 0)
                  {
                    ytop[s] = rtop;
                    yneg[s] = 0;
                  }
              }
            else
              {
                const int top = gmag_mul<WL, LF, GB>(Y[s], R, rtop, V[s], vtop[s], 0);
                if(g == 0)
                  {
                    ytop[s] = top;
                    yneg[s] = vneg[s];
                  }
              }
            // signed store: two's complement of a negative magnitude is ~v + 1,
            // i.e. ~v above the lowest nonzero limb, -v at it, 0 below
            const unsigned mask = group_mask<GB>();
            __syncwarp(mask);
            const bool ng = s != k && vneg[s];
            int low = WL;
            for(int w = g; w < WL; w += GB)
              if(Y[s][w])
                low = min(low, w);
            low = int(__reduce_min_sync(mask, unsigned(low)));
            u32 *dst = s == k ? at(XT, k, k) : s < k ? at(XT, k, s) : at(T, s, k);
            for(int w = g; w < WL; w += GB)
              {
                const u32 y = Y[s][w];
                dst[w] = !ng ? y : w < low ? 0u : w == low ? 0u - y : ~y;
              }
          }
        else if(tid < BT)
          {
            const int s = tid;
            u32 v[WL];
            if(s == k)
              {
#pragma unroll
                for(int w = 0; w < WL; ++w)
                  v[w] = R[w];
                store_signed<WL>(at(XT, k, k), v, false);
                ytop[s] = rtop;
                yneg[s] = 0;
              }
            else
              {
                mag_mul<WL, LF>(v, R, rtop, V[s], vtop[s]);
                const bool ng = vneg[s] != 0;
                store_signed<WL>(s < k ? at(XT, k, s) : at(T, s, k), v, ng);
                ytop[s] = top_reg<WL>(v);
                yneg[s] = ng;
              }
#pragma unroll
            for(int w = 0; w < WL; ++w)
              Y[s][w] = v[w];
          }
        __syncthreads();
        if(prof && tid == 0)
          {
            const long long c1 = clock64();
            cyc[1] += c1 - c0;
            c0 = c1;
          }
        // ---- C: T_is -= Y_i Y_s for i > k, 0 <= s <= i
        const int q0 = ((k + 1) * (k + 2)) / 2,
                  ntask = (BT * (BT + 1)) / 2 - q0;
        // with LA, warp 0 takes task (k+1, k+1) (index k+1) and the next
        // pivot; the other warps the rest
        const int p0 = LA ? tid - 32 : tid, pst = LA ? NT - 32 : NT,
                  nt2 = LA ? ntask - 1 : ntask;
        if(LA && tid < 32)
          {
            if(k + 1 < BT)
              next_pivot<WL, LF, BT>(at(T, k + 1, k + 1), Y[k + 1], ytop[k + 1],
                                     R, &rtop, Dm, T1, T2);
          }
        else
        for(int p = p0; p < nt2; p += pst)
          {
            const int q = (LA && p >= k + 1 ? p + 1 : p) + q0;
            int i = int((sqrtf(8.0f * q + 1.0f) - 1.0f) * 0.5f);
            while((i + 1) * (i + 2) / 2 <= q)
              ++i;
            while(i * (i + 1) / 2 > q)
              --i;
            const int s = q - i * (i + 1) / 2;
            u32 v[WL];
            // Y_i = l_ik is an entry of L (|.| <= 1): LF + 1 limbs; Y_s may be a
            // large entry of the inverse
            if(ytop[i] < LF + 1 && ytop[s] < LF + 1)
              mag_mul_small<WL, LF, LF + 1>(v, Y[i], Y[s]);
            else if(ytop[i] < LF + 1)
              mag_mul_small<WL, LF, LF + 1, WL>(v, Y[i], Y[s]);
            else
              mag_mul<WL, LF>(v, Y[i], ytop[i], Y[s], ytop[s]);
            rmw<WL>(s <= k ? at(XT, i, s) : at(T, i, s), v,
                    (yneg[i] ^ yneg[s]) != 0);
          }
        __syncthreads();
        if(prof && tid == 0)
          {
            const long long c1 = clock64();
            cyc[2] += c1 - c0;
          }
      }
    if(prof && tid == 0)
      for(int p = 0; p < 3; ++p)
        prof[blockIdx.x * 3 + p] = cyc[p];
  }
}
