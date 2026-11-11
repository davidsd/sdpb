// GPU backend for SDPB: batched block Cholesky with explicit inverses, and
// P = L^-1 B as products.  See Block_Factor.hxx.  Based on SDPviBe by Antoine
// Tilloy (https://github.com/AntoineTilloy/SDPviBe), the prototype of this
// backend: sdpvibe/linalg.py (chol_inv), sdpvibe/mpf.py (Chol: power-of-two
// equilibration, mpf_add), sdpvibe/basecase.py (leaves).  Portions adapted from
// SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).

#include "Block_Factor.hxx"
#include "basecase.cuh"
#include "leaf_fused.cuh"
#include "../rns/Limb_Sizes.hxx"
#include "../rns/Rns_Basis.hxx"
#include "../rns/rns_kernels.cuh"

#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <vector>

#define BF_CHECK(x)                                                           \
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
    constexpr int TH = 256;
    constexpr int LM = 64; // max limbs of a float
    constexpr int WM = 96; // max limbs of a wide integer
    int64_t nb(int64_t n)
    { return (n + TH - 1) / TH; }

    __device__ void negate(uint32_t *a, int n)
    {
      unsigned long long c = 1;
      for(int k = 0; k < n; ++k)
        {
          unsigned long long t = (unsigned long long)(~a[k]) + c;
          a[k] = uint32_t(t);
          c = t >> 32;
        }
    }
    __device__ int bitlen(const uint32_t *a, int n)
    {
      for(int k = n - 1; k >= 0; --k)
        if(a[k])
          return 32 * k + 32 - __clz(a[k]);
      return 0;
    }
    // dst[0..nd) = round_half_up(src / 2^s) for magnitudes (s > 0), or
    // src * 2^-s (s <= 0) (SDPviBe mpf.py shift_mag)
    __device__ void
    shift_mag(const uint32_t *src, int ns, long long s, uint32_t *dst, int nd)
    {
      if(s > 0)
        {
          long long q = s >> 5;
          int r = int(s & 31);
          for(int k = 0; k < nd; ++k)
            {
              long long i0 = q + k, i1 = q + k + 1;
              uint32_t lo = i0 < ns ? src[i0] : 0, hi = i1 < ns ? src[i1] : 0;
              dst[k] = r ? ((lo >> r) | (hi << (32 - r))) : lo;
            }
          long long sb = s - 1, qb = sb >> 5;
          int rb = int(sb & 31);
          unsigned long long c = (qb < ns) ? ((src[qb] >> rb) & 1u) : 0;
          for(int k = 0; k < nd && c; ++k)
            {
              unsigned long long t = (unsigned long long)dst[k] + c;
              dst[k] = uint32_t(t);
              c = t >> 32;
            }
        }
      else
        {
          long long t = -s, q = t >> 5;
          int r = int(t & 31);
          for(int k = 0; k < nd; ++k)
            {
              long long i0 = k - q, i1 = k - q - 1;
              uint32_t lo = (i0 >= 0 && i0 < ns) ? src[i0] : 0,
                       lo2 = (i1 >= 0 && i1 < ns) ? src[i1] : 0;
              dst[k] = r ? ((lo << r) | (lo2 >> (32 - r))) : lo;
            }
        }
    }
    // normalized float from a two's complement integer a (n limbs) * 2^base
    __device__ void
    store_norm(uint32_t *a, int n, long long base, int L, int32_t *m,
               int64_t plane, int64_t idx, int64_t *e)
    {
      bool ng = int32_t(a[n - 1]) < 0;
      if(ng)
        negate(a, n);
      int bl = bitlen(a, n);
      if(bl == 0)
        {
          for(int t = 0; t < L; ++t)
            m[(int64_t)t * plane + idx] = 0;
          e[idx] = MPF_EZ;
          return;
        }
      const int NBb = 32 * L - 2;
      long long s = bl - NBb;
      uint32_t r[LM + 1];
      shift_mag(a, n, s, r, L + 1);
      if(bitlen(r, L + 1) > NBb) // rounding carried to 2^NB
        {
          shift_mag(a, n, s + 1, r, L + 1);
          s += 1;
        }
      if(ng)
        negate(r, L);
      for(int t = 0; t < L; ++t)
        m[(int64_t)t * plane + idx] = int32_t(r[t]);
      e[idx] = base + s;
    }

    // one thread block per 16 x 16 tile, so that it reads and writes 16 lines
    // per limb plane instead of reading 16 and writing 256 (grid: x = tile
    // columns, y = tile rows, z = batch)
    __global__ void
    k_transpose(const int32_t *__restrict__ am, const int64_t *__restrict__ ae,
                int L, int batch, int64_t R, int64_t C, int64_t apl, int64_t abs,
                int64_t ald, int32_t *__restrict__ tm, int64_t *__restrict__ te)
    {
      const int64_t i = blockIdx.y * 16 + threadIdx.x % 16,
                    j = blockIdx.x * 16 + threadIdx.x / 16;
      if(i >= R || j >= C)
        return;
      const int64_t tot = int64_t(batch) * R * C, b = blockIdx.z;
      const int64_t idx = b * abs + j * ald + i; // A(i, j)
      const int64_t out = b * R * C + i * C + j; // T(j, i), T has C rows
      te[out] = ae[idx];
      for(int t = 0; t < L; ++t)
        tm[(int64_t)t * tot + out] = am[(int64_t)t * apl + idx];
    }

    __global__ void
    k_sub(const int32_t *__restrict__ sm, const int64_t *__restrict__ se,
          int64_t splane, int64_t R, int64_t C, int64_t r0, int64_t c0,
          int32_t *__restrict__ dm, int64_t *__restrict__ de, int64_t dplane,
          int64_t h, int64_t w, int L, int batch, bool insert)
    {
      const int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(idx >= int64_t(batch) * h * w)
        return;
      const int64_t b = idx / (h * w), rem = idx % (h * w);
      const int64_t ii = rem % h, jj = rem / h;
      const int64_t sidx = b * R * C + (c0 + jj) * R + r0 + ii;
      if(!insert)
        {
          de[idx] = se[sidx];
          for(int t = 0; t < L; ++t)
            dm[(int64_t)t * dplane + idx] = sm[(int64_t)t * splane + sidx];
        }
      else
        {
          // here (dm, de) is the small block and (sm, se) the large one
          const_cast<int64_t *>(se)[sidx] = de[idx];
          for(int t = 0; t < L; ++t)
            const_cast<int32_t *>(sm)[(int64_t)t * splane + sidx]
              = dm[(int64_t)t * dplane + idx];
        }
    }

    // C = A + sgn B, normalized (SDPviBe mpf.py mpf_add)
    // grid: x = columns, y = slices of 256 rows, z = blocks. Each operand has
    // its own layout. lower: entries above the diagonal are left alone.
    struct Layout
    {
      int64_t plane, bstride, ld;
    };
    __global__ void
    k_add(const int32_t *Am, const int64_t *Ae, Layout la, const int32_t *Bm,
          const int64_t *Be, Layout lb, int sgn, int L, int64_t rows, int32_t *Om,
          int64_t *Oe, Layout lo, bool lower)
    {
      const int64_t row = blockIdx.y * int64_t(blockDim.x) + threadIdx.x,
                    col = blockIdx.x, z = blockIdx.z;
      if(row >= rows || (lower && row < col))
        return;
      const int64_t ia = z * la.bstride + col * la.ld + row,
                    ib = z * lb.bstride + col * lb.ld + row,
                    io = z * lo.bstride + col * lo.ld + row;
      uint32_t a[LM], b[LM], wa[LM + 2], wb[LM + 2];
      for(int k = 0; k < L; ++k)
        {
          a[k] = uint32_t(Am[(int64_t)k * la.plane + ia]);
          b[k] = uint32_t(Bm[(int64_t)k * lb.plane + ib]);
        }
      bool na = a[L - 1] >> 31, nbg = ((b[L - 1] >> 31) != 0) ^ (sgn < 0);
      if(a[L - 1] >> 31)
        negate(a, L);
      if(b[L - 1] >> 31)
        negate(b, L);
      long long ea = Ae[ia], eb = Be[ib];
      if(bitlen(a, L) == 0)
        ea = MPF_EZ;
      if(bitlen(b, L) == 0)
        eb = MPF_EZ;
      if(ea == MPF_EZ && eb == MPF_EZ)
        {
          for(int k = 0; k < L; ++k)
            Om[(int64_t)k * lo.plane + io] = 0;
          Oe[io] = MPF_EZ;
          return;
        }
      long long t = (ea > eb ? ea : eb) - 32;
      long long da = ea - t, db = eb - t;
      long long sa = da < -(long long)(32 * (L + 3)) ? 32 * (L + 3) : -da;
      long long sb2 = db < -(long long)(32 * (L + 3)) ? 32 * (L + 3) : -db;
      if(ea == MPF_EZ)
        for(int k = 0; k < L + 2; ++k)
          wa[k] = 0;
      else
        shift_mag(a, L, sa, wa, L + 2);
      if(eb == MPF_EZ)
        for(int k = 0; k < L + 2; ++k)
          wb[k] = 0;
      else
        shift_mag(b, L, sb2, wb, L + 2);
      if(na)
        negate(wa, L + 2);
      if(nbg)
        negate(wb, L + 2);
      unsigned long long c = 0;
      for(int k = 0; k < L + 2; ++k)
        {
          unsigned long long v = (unsigned long long)wa[k] + wb[k] + c;
          wa[k] = uint32_t(v);
          c = v >> 32;
        }
      store_norm(wa, L + 2, t, L, Om, lo.plane, io, Oe);
    }

    __global__ void k_zero_e(int64_t *e, int64_t n)
    {
      const int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(i < n)
        e[i] = MPF_EZ;
    }

    __global__ void
    k_max_exp(const int64_t *e, int64_t rows, int64_t cols, int batch,
              int64_t bstride, int64_t ld, unsigned long long *out)
    {
      __shared__ long long v[256];
      long long mx = LLONG_MIN;
      const int64_t count = int64_t(batch) * rows * cols, per = rows * cols;
      for(int64_t q = int64_t(blockIdx.x) * blockDim.x + threadIdx.x;
          q < count; q += int64_t(gridDim.x) * blockDim.x)
        {
          const int64_t rem = q % per, idx = (q / per) * bstride
                                             + (rem / rows) * ld + rem % rows;
          if(e[idx] != MPF_EZ)
            mx = max(mx, (long long)e[idx]);
        }
      v[threadIdx.x] = mx;
      __syncthreads();
      for(int d = 128; d; d /= 2)
        {
          if(threadIdx.x < d)
            v[threadIdx.x] = max(v[threadIdx.x], v[threadIdx.x + d]);
          __syncthreads();
        }
      if(threadIdx.x == 0)
        atomicMax(out, (unsigned long long)v[0] ^ 0x8000000000000000ULL);
    }

    // flags[j*buckets + i/16] = 1 for every nonzero entry (i, j) of any block
    // (grid: x = columns, y = slices of TH rows, z = blocks)
    __global__ void k_column_rows(const int64_t *__restrict__ e, int64_t R,
                                  int64_t C, int64_t buckets,
                                  uint8_t *__restrict__ flags)
    {
      const int64_t i = int64_t(blockIdx.y) * blockDim.x + threadIdx.x;
      if(i >= R)
        return;
      const int64_t j = blockIdx.x;
      if(e[int64_t(blockIdx.z) * R * C + j * R + i] != MPF_EZ)
        flags[j * buckets + i / 16] = 1;
    }

    // diagnostic: per block, the largest exponent of its entries
    __global__ void k_block_max_exp(const int64_t *e, int64_t rows, int64_t cols,
                                    int batch, int64_t bstride, int64_t ld,
                                    unsigned long long *out)
    {
      const int64_t q = int64_t(blockIdx.x) * blockDim.x + threadIdx.x,
                    per = rows * cols;
      if(q >= batch * per)
        return;
      const int64_t rem = q % per, b = q / per,
                    idx = b * bstride + (rem / rows) * ld + rem % rows;
      if(e[idx] != MPF_EZ)
        atomicMax(out + b,
                  (unsigned long long)(long long)e[idx] ^ 0x8000000000000000ULL);
    }

    __global__ void k_unpack_symmetric(const int32_t *pm, const int64_t *pe,
                                       int32_t *m, int64_t *e, int n, int L)
    {
      int64_t q = int64_t(blockIdx.x) * blockDim.x + threadIdx.x;
      if(q >= int64_t(n) * n)
        return;
      int64_t i = q % n, j = q / n, hi = max(i, j), lo = min(i, j);
      int64_t p = lo * n - lo * (lo - 1) / 2 + hi - lo,
              packed = int64_t(n) * (n + 1) / 2;
      e[q] = pe[p];
      for(int t = 0; t < L; t++)
        m[int64_t(t) * n * n + q] = pm[int64_t(t) * packed + p];
    }

    // Exact power-of-two congruence. Mantissas are unchanged: only exponents
    // move. k_i = ceil(top(A_ii)/2), so (D A D)_ii < 1, D_ii = 2^-k_i.
    __global__ void k_equilibrate(const int64_t *ae, int64_t *be,
                                  int64_t *scale, int64_t n, int batch, int L)
    {
      int64_t q = int64_t(blockIdx.x) * blockDim.x + threadIdx.x;
      if(q >= int64_t(batch) * n * n)
        return;
      int64_t b = q / (n * n), i = q % n, j = (q / n) % n;
      int64_t ki = (ae[b * n * n + i * n + i] + 32LL * L - 1) >> 1;
      int64_t kj = (ae[b * n * n + j * n + j] + 32LL * L - 1) >> 1;
      if(j == 0)
        scale[b * n + i] = ki;
      be[q] = ae[q] == MPF_EZ ? MPF_EZ : ae[q] - ki - kj;
    }
    __global__ void k_unscale(int64_t *le, int64_t *xe, const int64_t *scale,
                              int64_t n, int batch)
    {
      int64_t q = int64_t(blockIdx.x) * blockDim.x + threadIdx.x;
      if(q >= int64_t(batch) * n * n)
        return;
      int64_t b = q / (n * n), i = q % n, j = (q / n) % n;
      // chol(A) = D^-1 chol(DAD); chol(A)^-1 = chol(DAD)^-1 D.
      if(le[q] != MPF_EZ)
        le[q] += scale[b * n + i];
      if(xe[q] != MPF_EZ)
        xe[q] -= scale[b * n + j];
    }

    // leaves: floats (batch tiles of 32 x 32) -> equilibrated wide fixed point
    // (row-major per tile, WL limbs, F = 32 LF fractional bits), row scales k_i
    // with a_ii < 2^(2 k_i)
    __global__ void
    k_leaf_in(const int32_t *__restrict__ m, const int64_t *__restrict__ e,
              int L, int tiles, int WL, int LF, int64_t pl, int64_t bs, int64_t ld,
              uint32_t *__restrict__ W, int64_t *__restrict__ kk)
    {
      const int64_t q = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(q >= int64_t(tiles) * 1024)
        return;
      const int64_t b = q / 1024, rem = q % 1024;
      const int i = int(rem % 32), j = int(rem / 32); // column-major input
      const int64_t idx = b * bs + j * ld + i;
      // the leaf reads the lower triangle only
      uint32_t *out = W + ((b * 32 + i) * 32 + j) * WL;
      if(i < j)
        {
          for(int t = 0; t < WL; ++t)
            out[t] = 0;
          return;
        }
      auto kdiag = [&](int r) -> long long {
        const int64_t d = b * bs + r * ld + r;
        // normalized: |m| has 32L-2 bits; a_rr < 2^(e + 32L - 2)
        long long top = e[d] + 32LL * L - 2;
        return (top + 1) >> 1; // ceil(top / 2), floor for negatives
      };
      const long long ki = kdiag(i), kj = kdiag(j);
      if(j == 0)
        kk[b * 32 + i] = ki;
      uint32_t a[LM], w[WM];
      for(int t = 0; t < L; ++t)
        a[t] = uint32_t(m[(int64_t)t * pl + idx]);
      if(e[idx] == MPF_EZ)
        {
          for(int t = 0; t < WL; ++t)
            out[t] = 0;
          return;
        }
      bool ng = int32_t(a[L - 1]) < 0;
      if(ng)
        negate(a, L);
      // value 2^(e - ki - kj) * a  ->  integer * 2^-F
      const long long s = -(e[idx] - ki - kj + 32LL * LF);
      shift_mag(a, L, s, w, WL);
      if(ng)
        negate(w, WL);
      for(int t = 0; t < WL; ++t)
        out[t] = w[t];
    }

    // leaf outputs: wide fixed point (row-major per tile) -> floats
    // (column-major), L_ij = l'_ij 2^(k_i) and X_ij = x'_ij 2^(-k_j), i >= j
    __global__ void
    k_leaf_out(const uint32_t *__restrict__ W, int WL, int LF,
               const int64_t *__restrict__ kk, int tiles, bool inverse, int L,
               int64_t pl, int64_t bs, int64_t ld, int32_t *__restrict__ m,
               int64_t *__restrict__ e)
    {
      const int64_t q = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(q >= int64_t(tiles) * 1024)
        return;
      const int64_t b = q / 1024, rem = q % 1024;
      const int i = int(rem % 32), j = int(rem / 32);
      const int64_t idx = b * bs + j * ld + i;
      uint32_t w[WM];
      if(i < j)
        {
          for(int t = 0; t < L; ++t)
            m[(int64_t)t * pl + idx] = 0;
          e[idx] = MPF_EZ;
          return;
        }
      const uint32_t *in = W + ((b * 32 + i) * 32 + j) * WL;
      for(int t = 0; t < WL; ++t)
        w[t] = in[t];
      const long long sc = inverse ? -kk[b * 32 + j] : kk[b * 32 + i];
      store_norm(w, WL, -32LL * LF + sc, L, m, pl, idx, e);
    }
  }

  namespace
  {
    template <int LF>
    void launch_leaf(uint32_t *W, uint32_t *Xw, uint32_t *R, uint32_t *S,
                     int tiles, const Gram_Options &o, cudaStream_t s)
    {
      constexpr int WL = leaf_limbs(LF);
      if(o.fused_leaf && o.leaf_lookahead)
        leaf_chol_inv<WL, LF, 32, 512, true, true>
          <<<tiles, 512, 0, s>>>(W, Xw, nullptr);
      else if(o.fused_leaf)
        leaf_chol_inv<WL, LF, 32, 512, true, false>
          <<<tiles, 512, 0, s>>>(W, Xw, nullptr);
      else
        {
          chol_tiles<WL, LF, 32><<<tiles, 256, 0, s>>>(W, R);
          BF_CHECK(cudaGetLastError());
          trinv_tiles<WL, LF, 32><<<tiles, 256, 0, s>>>(W, R, Xw, S);
        }
      BF_CHECK(cudaGetLastError());
    }
  }

  // ---------------------------------------------------------------- helpers
  void mpf_transpose(const Mpf_Batch &A, Mpf_Batch &T, cudaStream_t s)
  {
    k_transpose<<<dim3(unsigned((A.cols + 15) / 16), unsigned((A.rows + 15) / 16),
                       unsigned(A.batch)),
                  TH, 0, s>>>(A.m, A.e, A.L, A.batch, A.rows, A.cols, A.pl(),
                              A.bs(), A.lda(), T.m, T.e);
    BF_CHECK(cudaGetLastError());
  }
  void mpf_sub_block(const Mpf_Batch &src, int64_t r0, int64_t c0,
                     Mpf_Batch &dst, bool insert, cudaStream_t s)
  {
    k_sub<<<nb(dst.entries()), TH, 0, s>>>(
      src.m, src.e, src.entries(), src.rows, src.cols, r0, c0, dst.m, dst.e,
      dst.entries(), dst.rows, dst.cols, src.L, src.batch, insert);
    BF_CHECK(cudaGetLastError());
  }
  void mpf_batch_add(const Mpf_Batch &A, const Mpf_Batch &B, int sgn,
                     const Mpf_Batch &C, cudaStream_t s, bool lower_only)
  {
    k_add<<<dim3(unsigned(A.cols), unsigned((A.rows + TH - 1) / TH),
                 unsigned(A.batch)),
            TH, 0, s>>>(A.m, A.e, Layout{A.pl(), A.bs(), A.lda()}, B.m, B.e,
                        Layout{B.pl(), B.bs(), B.lda()}, sgn, A.L, A.rows, C.m,
                        C.e, Layout{C.pl(), C.bs(), C.lda()}, lower_only);
    BF_CHECK(cudaGetLastError());
  }
  void mpf_unpack_symmetric(const Mpf_Batch &packed, Mpf_Batch &M,
                            cudaStream_t s)
  {
    k_unpack_symmetric<<<nb(M.entries()), TH, 0, s>>>(packed.m, packed.e, M.m,
                                                      M.e, int(M.rows), M.L);
    BF_CHECK(cudaGetLastError());
  }
  void mpf_zero(Mpf_Batch &A, cudaStream_t s)
  {
    BF_CHECK(cudaMemsetAsync(A.m, 0, size_t(A.L) * A.entries() * 4, s));
    k_zero_e<<<nb(A.entries()), TH, 0, s>>>(A.e, A.entries());
    BF_CHECK(cudaGetLastError());
  }

  // ---------------------------------------------------------------- class
  Block_Factor::Block_Factor(int n, int batch, int L, int digits,
                             Gram_Options options, cudaStream_t stream)
      : n_(n), batch_(batch), L_(L), D_(digits), opt_(options), stream_(stream)
  {
    if(n % 32 != 0)
      throw std::invalid_argument("Block_Factor: n % 32 != 0 (pad)");
    if(rounded_L(L) == 0)
      throw std::invalid_argument(
        "Block_Factor: L above the compiled sizes (rns/Limb_Sizes.hxx)");
    // sizes for the planner in gemm_for. The basis is the largest a product of
    // this shape can use; these margins give 132, 144, 264 and 60 at L = 28.
    {
      const Rns_Basis basis = Rns_Basis::for_products(D_, int64_t(1) << 16);
      plan_primes_ = double(basis.primes.size()) + 2;
      plan_padded_ = basis.num_padded;
      plan_crt_digits_ = basis.crt_digits_128();
      plan_limbs_ = basis.crt_digits / 4.0 + 2;
    }
    // options.workspace: the caller's, shared with other engines on this stream
    if(opt_.share_workspace && !opt_.workspace)
      {
        workspace_ = std::make_unique<Rns_Workspace>(stream_);
        opt_.workspace = workspace_.get();
      }
    BF_CHECK(cudaEventCreate(&ev0_));
    BF_CHECK(cudaEventCreate(&ev1_));
  }

  Block_Factor::~Block_Factor()
  {
    pool_flush();
    for(void *p : {(void *)leaf_W_, (void *)leaf_Xw_, (void *)leaf_S_,
                   (void *)leaf_R_, (void *)leaf_kk_})
      if(p)
        cudaFree(p);
    gemms_.clear();
    grams_.clear();
    cudaEventDestroy(ev0_);
    cudaEventDestroy(ev1_);
  }

  void Block_Factor::tic()
  {
    if(opt_.timing)
      BF_CHECK(cudaEventRecord(ev0_, stream_));
  }
  void Block_Factor::toc(const std::string &phase)
  {
    if(!opt_.timing)
      return;
    BF_CHECK(cudaEventRecord(ev1_, stream_));
    BF_CHECK(cudaEventSynchronize(ev1_));
    float ms = 0;
    BF_CHECK(cudaEventElapsedTime(&ms, ev0_, ev1_));
    ms_[phase] += ms;
  }

  namespace
  {
    // thrown by refined_panel, caught by factor(): the estimate passed the limit
    struct Correction_Needed
    {
    };

    // host time of a scope, added to *ms when ms != nullptr
    struct Host_Timer
    {
      double *ms;
      std::chrono::steady_clock::time_point t0;
      explicit Host_Timer(double *p)
          : ms(p), t0(std::chrono::steady_clock::now()) {}
      ~Host_Timer()
      {
        if(ms)
          *ms += std::chrono::duration<double, std::milli>(
                   std::chrono::steady_clock::now() - t0)
                   .count();
      }
    };
  }

  Mpf_Batch Block_Factor::alloc(int64_t rows, int64_t cols)
  {
    Host_Timer ht(opt_.timing ? &ms_["host alloc and free"] : nullptr);
    Mpf_Batch M;
    M.L = L_;
    M.batch = batch_;
    M.rows = rows;
    M.cols = cols;
    auto allocate = [&](void **p, size_t bytes) {
      auto error = malloc_async_retry(p, bytes, stream_);
      if(error == cudaErrorMemoryAllocation && !pool_.empty())
        {
          (void)cudaGetLastError();
          pool_flush();
          error = malloc_async_retry(p, bytes, stream_);
        }
      if(error == cudaErrorMemoryAllocation)
        {
          // No cached plan is borrowed across alloc(). Wait for all queued
          // products before releasing their buffers, then retry the temporary
          // allocation. Product creation already has a separate memory budget.
          (void)cudaGetLastError();
          check();
        }
      if(error == cudaErrorMemoryAllocation && opt_.workspace)
        {
          // products on the stream are done (check() above); they rebuild the
          // workspace when next used
          opt_.workspace->release();
          ++allocation_evictions_;
          error = malloc_async_retry(p, bytes, stream_);
        }
      while(error == cudaErrorMemoryAllocation
            && (!gemms_.empty() || !grams_.empty()))
        {
          if(!gemms_.empty())
            {
              auto big = gemms_.begin();
              for(auto it = gemms_.begin(); it != gemms_.end(); ++it)
                if(bytes_[it->first] > bytes_[big->first])
                  big = it;
              bytes_.erase(big->first);
              gemms_.erase(big);
            }
          else
            grams_.clear();
          ++allocation_evictions_;
          error = malloc_async_retry(p, bytes, stream_);
        }
      return error;
    };
    const size_t bm = size_t(L_) * M.entries() * 4, be = size_t(M.entries()) * 8;
    M.m = static_cast<int32_t *>(pool_get(bm));
    M.e = static_cast<int64_t *>(pool_get(be));
    if(!M.m)
      BF_CHECK(allocate(reinterpret_cast<void **>(&M.m), bm));
    if(!M.e)
      {
        auto error = allocate(reinterpret_cast<void **>(&M.e), be);
        if(error != cudaSuccess)
          {
            BF_CHECK(cudaFreeAsync(M.m, stream_));
            BF_CHECK(error);
          }
      }
    if(pool_active_)
      {
        live_[M.m] = bm;
        live_[M.e] = be;
      }
    return M;
  }

  void *Block_Factor::pool_get(size_t bytes)
  {
    if(!pool_active_)
      return nullptr;
    auto it = pool_.find(bytes);
    if(it == pool_.end())
      return nullptr;
    void *p = it->second;
    pool_.erase(it);
    return p;
  }

  void Block_Factor::pool_put(void *p, size_t bytes)
  {
    pool_.emplace(bytes, p);
  }

  void Block_Factor::pool_flush()
  {
    for(auto &kv : pool_)
      cudaFreeAsync(kv.second, stream_);
    pool_.clear();
  }
  void Block_Factor::release(Mpf_Batch &M)
  {
    Host_Timer ht(opt_.timing ? &ms_["host alloc and free"] : nullptr);
    if(pool_active_ && M.m)
      {
        live_.erase(M.m);
        live_.erase(M.e);
        pool_put(M.m, size_t(M.L) * M.entries() * 4);
        pool_put(M.e, size_t(M.entries()) * 8);
      }
    else
      {
        if(M.m)
          BF_CHECK(cudaFreeAsync(M.m, stream_));
        if(M.e)
          BF_CHECK(cudaFreeAsync(M.e, stream_));
      }
    M.m = nullptr;
    M.e = nullptr;
  }

  Mpf_Gemm &Block_Factor::gemm_for(int m, int n, int k, int digits)
  {
    auto key = std::make_tuple(m, n, k, digits, product_guard_);
    auto it = gemms_.find(key);
    if(it == gemms_.end())
      {
        Gram_Options o = opt_;
        if(product_guard_ > 0)
          o.crt_guard_bits = product_guard_;
        o.timing = opt_.profile_products || opt_.product_timing;
        // no host synchronization per product: checked at the end of
        // factor() and solve()
        o.defer_crt_check = true;
        // small products are latency-bound: all primes in one group (one
        // residue launch per operand, one batched GEMM, one reduction)
        if(double(batch_) * k * (m + std::min(n, 1024)) <= double(1 << 21))
          {
            o.prime_group = 256;
            o.max_entries_per_residue_batch = int64_t(1) << 16;
          }
        // device memory: choose the panel width and whether to cache the
        // first operand's residues so that the product's buffers (as
        // Rns_Gemm and Mpf_Gemm allocate them) take at most 60% of the free
        // memory; if they don't fit, free the other cached products first
        // (synchronize first: memory freed with cudaFreeAsync goes back to the
        // device only at a synchronization, and Rns_Gemm uses cudaMalloc)
        pool_flush();
        BF_CHECK(cudaStreamSynchronize(stream_));
        // A caller may retain an async allocation pool. Its unused reserve is
        // not available to cudaMalloc and otherwise looks like live memory to
        // the planner, causing avoidable eviction of reusable product plans.
        int device;
        cudaMemPool_t pool;
        BF_CHECK(cudaGetDevice(&device));
        BF_CHECK(cudaDeviceGetDefaultMemPool(&pool, device));
        BF_CHECK(cudaMemPoolTrimTo(pool, 0));
        const double r = plan_primes_, rp = plan_padded_, g = o.prime_group,
                     bm = double(batch_) * m;
        // residues of C: uint16 per prime, or with reduce_to_crt (groups of
        // at most 32 primes) the CRT rows (2 rp bytes) and a float
        const bool fused
          = o.reduce_to_crt && !o.fused_crt && o.prime_group <= 32;
        auto need = [&](int pn, bool ca) {
          const double crt = std::min(double(o.crt_batch), bm * pn);
          return 3.0 * g * batch_ * double(k) * pn   // RB
                 + 12.0 * g * bm * pn                  // int32 products
                 + (fused ? 2.0 * rp + 4 : 2.0 * r) * bm * pn // residues of C
                 + crt * ((fused ? 0 : 2 * rp) + 4 * plan_crt_digits_) // CRT buffers
                 + 4.0 * plan_limbs_ * bm * pn                 // exact panel
                 + (ca ? 3.0 * r * bm * k : 3.0 * g * bm * k); // residues of A
        };
        bool cache_a = o.cache_a;
        // SDPB_GPU_PANEL: the widest column panel (default 512: no bit changes
        // on five real dumps and 1.2-18 GB saved, at most 1.7% time)
        const char *pf = getenv("SDPB_GPU_PANEL");
        const int widest = pf ? atoi(pf) : 512;
        int panel = std::min(n, widest);
        for(int pass = 0; pass < 2; ++pass)
          {
            size_t free_b = 0, total_b = 0;
            BF_CHECK(cudaMemGetInfo(&free_b, &total_b));
            // fraction of the free memory one product may plan for
            // (SDPB_GPU_BUDGET, default 0.7: on the real nmax14 blocks 0.6
            // left the residue cache off and cost 1.6x in P = X B; with all 36
            // nmax18 blocks at once 0.8 ran out of memory)
            const char *bf = getenv("SDPB_GPU_BUDGET");
            const double budget
              = (bf ? atof(bf) : 0.7) * double(free_b)
                + (opt_.workspace ? double(opt_.workspace->bytes()) : 0.0);
            cache_a = o.cache_a;
            panel = std::min(n, widest);
            while(need(panel, cache_a) > budget && panel > 128)
              panel /= 2;
            if(need(panel, cache_a) > budget && cache_a)
              {
                cache_a = false;
                panel = std::min(n, widest);
                while(need(panel, cache_a) > budget && panel > 128)
                  panel /= 2;
              }
            if(getenv("SDPB_GPU_VERBOSE"))
              printf("  [Block_Factor] product m=%d n=%d k=%d batch=%d: panel "
                     "%d, cache_a %d, "
                     "buffers %.2f GB, free %.2f GB, %zu cached products\n",
                     m, n, k, batch_, panel, int(cache_a),
                     need(panel, cache_a) / 1e9, free_b / 1e9, gemms_.size());
            if(need(panel, cache_a) <= budget || gemms_.empty() || opt_.workspace)
              break;
            // free the largest cached products until this one fits
            check();
            double freed = 0;
            while(!gemms_.empty() && need(panel, cache_a) > budget + freed)
              {
                auto big = gemms_.begin();
                for(auto q = gemms_.begin(); q != gemms_.end(); ++q)
                  if(bytes_[q->first] > bytes_[big->first])
                    big = q;
                freed += bytes_[big->first];
                bytes_.erase(big->first);
                gemms_.erase(big);
              }
            BF_CHECK(cudaStreamSynchronize(stream_));
          }
        bytes_[key] = need(panel, cache_a);
        o.cache_a = cache_a;
        it = gemms_
               .emplace(key, std::make_unique<Mpf_Gemm>(
                               m, n, k, batch_, L_, digits, o, panel, stream_))
               .first;
      }
    return *it->second;
  }

  // C = op(U) op(V) via Mpf_Gemm's A'^T B' with A' = op(U)^T, B' = op(V)
  void Block_Factor::mm(const Mpf_Batch &U, bool tU, const Mpf_Batch &V,
                        bool tV, Mpf_Batch &C, int digits, int triangles)
  {
    if(digits == 0)
      digits = D_;
    const int64_t m = tU ? U.cols : U.rows, n = tV ? V.rows : V.cols;
    Mpf_Batch Ap = U, Bp = V;
    bool freeA = false, freeB = false;
    const bool same = triangles & Product_Same_Operands;
    if(same && (tU || !tV || U.m != V.m))
      throw std::invalid_argument("Block_Factor::mm: same operands are U U^T");
    if(!tU)
      {
        Ap = transposed(U);
        freeA = true;
      }
    if(same)
      Bp = Ap;
    else if(tV)
      {
        Bp = transposed(V);
        freeB = true;
      }
    if(!C.m)
      C = alloc(m, n);
    product(Ap, Bp, C, digits, triangles);
    if(freeA)
      release(Ap);
    if(freeB)
      release(Bp);
  }

  Mpf_Batch Block_Factor::transposed(const Mpf_Batch &M)
  {
    Mpf_Batch T = alloc(M.cols, M.rows);
    tic();
    mpf_transpose(M, T, stream_);
    toc("transpose");
    return T;
  }

  void Block_Factor::product(const Mpf_Batch &Ap, const Mpf_Batch &Bp,
                             const Mpf_Batch &C, int digits, int triangles)
  {
    if(digits == 0)
      digits = D_;
    const int64_t m = Ap.cols, k = Ap.rows, n = Bp.cols;
    tic();
    auto &g = gemm_for(int(m), int(n), int(k), digits);
    std::map<std::string, double> before, before2;
    if(opt_.profile_products)
      {
        before = g.timings();
        before2 = g.gemm().timings();
      }
    const int structure
      = triangles & (Product_Same_Operands | Product_Lower_Only | Product_Negate);
    if(C.ld)
      g.set_output_layout(C.pl(), C.bs(), C.lda());
    g.multiply(Ap.view(), Bp.view(), C.m, C.e,
               (skip_triangular_zeros_ ? triangles & 7 : 0) | structure);
    const std::string shape = "k=" + std::to_string(k) + " m=" + std::to_string(m)
                              + " n=" + std::to_string(n)
                              + " D=" + std::to_string(digits);
    toc("role/" + role_ + "/" + shape);
    if(opt_.timing)
      ms_["count/" + role_ + "/" + shape] += 1;
    if(opt_.profile_products)
      {
        for(const auto &kv : g.timings())
          ms_["kernel/" + role_ + "/" + kv.first] += kv.second - before[kv.first];
        for(const auto &kv : g.gemm().timings())
          ms_["kernel/" + role_ + "/" + kv.first] += kv.second - before2[kv.first];
      }
  }

  void Block_Factor::mm_flat(const Mpf_Batch &U, const Mpf_Batch &V, int32_t *d_m,
                             int64_t *d_e)
  {
    tic();
    gemm_for(int(U.cols), int(V.cols), int(U.rows), D_)
      .multiply(U.view(), V.view(), d_m, d_e, 0, true);
    toc(std::string("flat product"));
  }

  void Block_Factor::solve_flat(const Mpf_Batch &Xt, const Mpf_Batch &B,
                                int32_t *d_m, int64_t *d_e, int guard)
  {
    role_ = "solve P=X B";
    product_guard_ = guard;
    structured_product(Xt, B, d_m, d_e, true);
    product_guard_ = 0;
  }

  void Block_Factor::collect_product_timings(
    std::map<std::string, double> &ms) const
  {
    for(const auto &g : gemms_)
      {
        for(const auto &kv : g.second->timings())
          ms[kv.first] += kv.second;
        for(const auto &kv : g.second->gemm().timings())
          ms[kv.first] += kv.second;
      }
  }

  void Block_Factor::leaf(const Mpf_Batch &A, const Mpf_Batch &Lf,
                          const Mpf_Batch &X)
  {
    // the compiled size at or above L: more fractional bits than asked for if
    // L is not one of them (see rns/Limb_Sizes.hxx)
    const int LF = rounded_L(L_), WL = leaf_limbs(LF);
    const int tiles = A.batch;
    const size_t wb = size_t(tiles) * 1024 * WL * 4;
    if(!leaf_W_)
      {
        BF_CHECK(cudaStreamSynchronize(stream_));
        BF_CHECK(cudaMalloc(&leaf_W_, wb));
        BF_CHECK(cudaMalloc(&leaf_Xw_, wb));
        if(!opt_.fused_leaf)
          {
            BF_CHECK(cudaMalloc(&leaf_S_, wb));
            BF_CHECK(cudaMalloc(&leaf_R_, size_t(tiles) * 32 * WL * 4));
          }
        BF_CHECK(cudaMalloc(&leaf_kk_, size_t(tiles) * 32 * 8));
      }
    uint32_t *W = leaf_W_, *R = leaf_R_, *Xw = leaf_Xw_, *S = leaf_S_;
    int64_t *kk = leaf_kk_;
    tic();
    k_leaf_in<<<nb(int64_t(tiles) * 1024), TH, 0, stream_>>>(
      A.m, A.e, L_, tiles, WL, LF, A.pl(), A.bs(), A.lda(), W, kk);
    BF_CHECK(cudaGetLastError());
    toc("leaf: floats to fixed point");
    tic();
#define BF_LEAF(LT, ...)                                                      \
  case LT:                                                                    \
    launch_leaf<LT>(W, Xw, R, S, tiles, opt_, stream_);                       \
    break;
    switch(LF)
      {
        SDPB_GPU_EACH_L(BF_LEAF)
      default:
        throw std::logic_error("Block_Factor::leaf: not a compiled size");
      }
#undef BF_LEAF
    toc("leaf: Cholesky and inverse");
    tic();
    k_leaf_out<<<nb(int64_t(tiles) * 1024), TH, 0, stream_>>>(
      W, WL, LF, kk, tiles, false, L_, Lf.pl(), Lf.bs(), Lf.lda(), Lf.m, Lf.e);
    k_leaf_out<<<nb(int64_t(tiles) * 1024), TH, 0, stream_>>>(
      Xw, WL, LF, kk, tiles, true, L_, X.pl(), X.bs(), X.lda(), X.m, X.e);
    BF_CHECK(cudaGetLastError());
    toc("leaf: fixed point to floats");
    if(corrupt_leaf_limbs_ > 0)
      {
        // block 0 of the tile, limb planes t < corrupt_leaf_limbs_
        BF_CHECK(cudaStreamSynchronize(stream_));
        const std::vector<int32_t> none(32, 0);
        for(int t = 0; t < corrupt_leaf_limbs_; ++t)
          for(int j = 0; j < 32; ++j)
            BF_CHECK(cudaMemcpy(X.m + int64_t(t) * X.pl() + j * X.lda(), none.data(),
                                32 * 4, cudaMemcpyHostToDevice));
        printf("CORRUPTED one leaf: low %d of %d limbs of X, block 0\n",
               corrupt_leaf_limbs_, L_);
        corrupt_leaf_limbs_ = 0;
      }
  }

  void Block_Factor::refined_panel(const Mpf_Batch &B, const Mpf_Batch &Lf,
                                   const Mpf_Batch &X, const Mpf_Batch &Z,
                                   int refine, int cdigits)
  {
    role_ = "panel Z=A21 X11^T";
    // X11^T is used twice when the panel is corrected
    Mpf_Batch Xt = transposed(X), Bt = transposed(B);
    product(Bt, Xt, Z, 0, 4);
    release(Bt);
    if(log_panels_)
      {
        log_panel("A21", B);
        log_panel("X11", X);
        log_panel("Z", Z);
      }
    if(X.cols < refine_min_ || (refine_max_ && X.cols > refine_max_))
      {
        if(correction_limit_ == 0)
          refine = 0;
      }
    if(refine && correction_limit_ < 0)
      {
        tic();
        long long ta, tb;
        max_tops(B, X, ta, tb);
        int bits_k = 0;
        for(int64_t k = X.cols - 1; k > 0; k >>= 1)
          ++bits_k;
        const long long estimate
          = (ta == LLONG_MIN || tb == LLONG_MIN)
              ? LLONG_MIN
              : bits_k + 6 + ta + tb - 8LL * D_;
        toc("correction: estimate");
        if(log_panels_)
          printf("PANELDEC n1=%ld n2=%ld estimate %lld limit %d: %s\n",
                 long(X.cols), long(B.rows), estimate, correction_limit_,
                 estimate > correction_limit_ ? "restart" : "skip");
        if(estimate > correction_limit_)
          {
            release(Xt);
            throw Correction_Needed{};
          }
        if(X.cols < refine_min_ || (refine_max_ && X.cols > refine_max_)
           || (!refine_min_ && !refine_max_))
          {
            refine = 0;
            ++panels_skipped_;
          }
      }
    if(!refine && measure_residual_)
      {
        Mpf_Batch ZL, R = alloc(B.rows, B.cols);
        role_ = "diagnostic";
        mm(Z, false, Lf, true, ZL, 0, 4);
        mpf_batch_add(B, ZL, -1, R, stream_);
        log_panel("R", R);
        release(ZL);
        release(R);
      }
    for(int it = 0; it < refine; ++it)
      {
        Mpf_Batch ZL, R = alloc(B.rows, B.cols), C;
        role_ = "correction ZL=Z L11^T";
        mm(Z, false, Lf, true, ZL, 0, 4);
        tic();
        mpf_batch_add(B, ZL, -1, R, stream_);
        toc("correction: add");
        tic();
        const int digits = cdigits < 0 ? correction_digits(R, X) : cdigits;
        toc("correction: digits");
        role_ = "correction C=R X11^T";
        Mpf_Batch Rt = transposed(R);
        C = alloc(B.rows, B.cols);
        product(Rt, Xt, C, digits, 4);
        release(Rt);
        tic();
        mpf_batch_add(Z, C, 1, Z, stream_);
        toc("correction: add");
        if(log_panels_)
          log_panel("R", R);
        release(ZL);
        release(R);
        release(C);
      }
    release(Xt);
  }

  // diagnostic: log2 of the largest entry of M, per block (exponent + 32 L - 2)
  void Block_Factor::log_panel(const char *what, const Mpf_Batch &M)
  {
    unsigned long long *d;
    std::vector<unsigned long long> h(M.batch);
    BF_CHECK(cudaMallocAsync(&d, M.batch * 8, stream_));
    BF_CHECK(cudaMemsetAsync(d, 0, M.batch * 8, stream_));
    k_block_max_exp<<<nb(M.entries()), TH, 0, stream_>>>(
      M.e, M.rows, M.cols, M.batch, M.bs(), M.lda(), d);
    BF_CHECK(cudaGetLastError());
    BF_CHECK(cudaMemcpyAsync(h.data(), d, M.batch * 8, cudaMemcpyDeviceToHost,
                             stream_));
    BF_CHECK(cudaStreamSynchronize(stream_));
    BF_CHECK(cudaFreeAsync(d, stream_));
    printf("PANEL %s rows=%ld cols=%ld :", what, long(M.rows), long(M.cols));
    for(unsigned long long v : h)
      {
        const long long t = (long long)(v ^ 0x8000000000000000ULL);
        printf(" %lld", v == 0 ? -1000000LL : t + 32 * L_ - 2);
      }
    printf("\n");
  }

  void Block_Factor::symmetric_product(const Mpf_Batch &A, Mpf_Batch &C)
  {
    const int n = int(A.rows), k = int(A.cols);
    // alloc() may evict the Gram cache on memory pressure. Allocate every
    // temporary before borrowing a plan so eviction cannot invalidate it.
    auto At = alloc(k, n);
    const int64_t packed = int64_t(n) * (n + 1) / 2;
    auto P = alloc(packed, 1);
    if(!C.m)
      C = alloc(n, n);
    auto key = std::make_pair(n, k);
    auto it = grams_.find(key);
    if(it == grams_.end())
      {
        BF_CHECK(cudaStreamSynchronize(stream_));
        int device;
        cudaMemPool_t pool;
        BF_CHECK(cudaGetDevice(&device));
        BF_CHECK(cudaDeviceGetDefaultMemPool(&pool, device));
        BF_CHECK(cudaMemPoolTrimTo(pool, 0));
        Gram_Options o = opt_;
        o.timing = false;
        it = grams_
               .emplace(key,
                        std::make_unique<Mpf_Gram>(n, L_, D_, k, o, stream_))
               .first;
      }
    mpf_transpose(A, At, stream_);
    it->second->scan(At.view());
    it->second->add(At.view());
    it->second->finish(P.m, P.e);
    k_unpack_symmetric<<<nb(int64_t(n) * n), TH, 0, stream_>>>(P.m, P.e, C.m,
                                                               C.e, n, L_);
    BF_CHECK(cudaGetLastError());
    release(P);
    release(At);
  }

  // (Lf, X) of A (n x n): SDPviBe's chol_inv recursion, in place in Lf and X
  // (zero above the diagonal on entry); A's parts are read in place too
  void Block_Factor::rec(const Mpf_Batch &A, const Mpf_Batch &Lf,
                         const Mpf_Batch &X, int refine, int min_n, int cdigits,
                         int gram_min)
  {
    const int64_t n = A.rows;
    if(n == 32)
      {
        leaf(A, Lf, X);
        return;
      }
    const int64_t n1 = (n / 2 + 31) / 32 * 32, n2 = n - n1;
    const Mpf_Batch A21 = A.block(n1, 0, n2, n1), A22 = A.block(n1, n1, n2, n2);
    const Mpf_Batch L11 = Lf.block(0, 0, n1, n1), L21 = Lf.block(n1, 0, n2, n1),
                    L22 = Lf.block(n1, n1, n2, n2);
    const Mpf_Batch X11 = X.block(0, 0, n1, n1), X21 = X.block(n1, 0, n2, n1),
                    X22 = X.block(n1, n1, n2, n2);
    rec(A.block(0, 0, n1, n1), L11, X11, refine, min_n, cdigits, gram_min);
    refined_panel(A21, L11, X11, L21, refine, cdigits);
    role_ = "update W=L21 L21^T";
    Mpf_Batch L21t = transposed(L21); // for W and T
    Mpf_Batch Wm, T;
    tic();
    if(gram_min && n2 >= gram_min)
      {
        symmetric_product(L21, Wm);
        toc("role/" + role_ + "/packed");
      }
    else
      {
        Wm = alloc(n2, n2);
        product(L21t, L21t, Wm, 0,
                lower_update_ ? Product_Same_Operands | Product_Lower_Only : 0); // L21 L21^T
      }
    Mpf_Batch Sm = alloc(n2, n2);
    tic();
    // S = A22 - L21 L21^T. With lower_update, entries above the diagonal of
    // Sm are not written, and nothing below reads them.
    mpf_batch_add(A22, Wm, -1, Sm, stream_, lower_update_);
    toc("add: Schur complement");
    release(Wm);
    // L21 X11, straight after L21 L21^T: the same residues of L21^T
    role_ = "inverse T=L21 X11";
    T = alloc(n2, n1);
    if(reuse_residues_ && !(gram_min && n2 >= gram_min))
      {
        auto w = gemms_.find(std::make_tuple(int(n2), int(n2), int(n1), D_, product_guard_));
        if(w != gemms_.end())
          {
            auto &t = gemm_for(int(n2), int(n1), int(n1), D_);
            // gemm_for may have evicted the update's product
            w = gemms_.find(std::make_tuple(int(n2), int(n2), int(n1), D_, product_guard_));
            if(w != gemms_.end() && w->second.get() != &t)
              t.gemm().borrow_a_residues(w->second->gemm());
          }
      }
    product(L21t, X11, T, 0, 2);
    release(L21t);
    rec(Sm, L22, X22, refine, min_n, cdigits, gram_min);
    release(Sm);
    role_ = "inverse X21=X22 T";
    // X21 = -X22 L21 X11 (negated in the rounding)
    Mpf_Batch X22t = transposed(X22);
    product(X22t, T, X21, 0, 1 | Product_Negate);
    release(X22t);
    release(T);
  }

  // Congruence-scaled inputs have diagonal in [1/4,1), hence every exact
  // Cholesky row has norm < 1. Budget the absolute error of the SMALL correction
  // below 2^-(32 L + 4), 2^-900 at L = 28: k*64*2^(-8D)*max|R|*max|X| <= that.
  // The full-precision residual is untouched. Cap D at the available precision;
  // reaching that cap does not certify arbitrary ill-conditioned inputs.
  void Block_Factor::max_tops(const Mpf_Batch &A, const Mpf_Batch &B,
                              long long &ta, long long &tb)
  {
    unsigned long long *d, h[2] = {0, 0};
    BF_CHECK(cudaMallocAsync(&d, 16, stream_));
    BF_CHECK(cudaMemsetAsync(d, 0, 16, stream_));
    k_max_exp<<<256, 256, 0, stream_>>>(A.e, A.rows, A.cols, A.batch, A.bs(),
                                        A.lda(), d);
    k_max_exp<<<256, 256, 0, stream_>>>(B.e, B.rows, B.cols, B.batch, B.bs(),
                                        B.lda(), d + 1);
    BF_CHECK(cudaGetLastError());
    BF_CHECK(cudaMemcpyAsync(h, d, 16, cudaMemcpyDeviceToHost, stream_));
    BF_CHECK(cudaStreamSynchronize(stream_));
    BF_CHECK(cudaFreeAsync(d, stream_));
    const auto top = [&](unsigned long long v) {
      const long long e = (long long)(v ^ 0x8000000000000000ULL);
      return e == LLONG_MIN ? e : e + 32LL * L_ - 2;
    };
    ta = top(h[0]);
    tb = top(h[1]);
  }

  int Block_Factor::correction_digits(const Mpf_Batch &R, const Mpf_Batch &X)
  {
    unsigned long long *d, h[2] = {0, 0};
    BF_CHECK(cudaMallocAsync(&d, 16, stream_));
    BF_CHECK(cudaMemsetAsync(d, 0, 16, stream_));
    k_max_exp<<<256, 256, 0, stream_>>>(R.e, R.rows, R.cols, R.batch, R.bs(),
                                        R.lda(), d);
    k_max_exp<<<256, 256, 0, stream_>>>(X.e, X.rows, X.cols, X.batch, X.bs(),
                                        X.lda(), d + 1);
    BF_CHECK(cudaGetLastError());
    BF_CHECK(cudaMemcpyAsync(h, d, 16, cudaMemcpyDeviceToHost, stream_));
    BF_CHECK(cudaStreamSynchronize(stream_));
    BF_CHECK(cudaFreeAsync(d, stream_));
    long long er = (long long)(h[0] ^ 0x8000000000000000ULL);
    long long ex = (long long)(h[1] ^ 0x8000000000000000ULL);
    int bits_k = 0;
    for(int64_t k = R.cols - 1; k > 0; k >>= 1)
      ++bits_k;
    long long bits = er == LLONG_MIN || ex == LLONG_MIN
                       ? 0
                       : 32 * L_ + 4 + 6 + bits_k + er + ex + 2 * (32 * L_ - 2);
    // The INT8 layout requires a multiple of four base-256 digits. Round up
    // (never down), preserving the error budget.
    int digits
      = int(std::max(8LL, std::min((long long)D_, 4 * ((bits + 31) / 32))));
    if(getenv("SDPB_GPU_CORRECTION_VERBOSE"))
      printf("  correction %ld x %ld: digits %d, max exponents %lld %lld\n",
             long(R.rows), long(R.cols), digits, er + 32 * L_ - 2,
             ex + 32 * L_ - 2);
    return digits;
  }

  void Block_Factor::check()
  {
    for(auto &g : gemms_)
      g.second->gemm().check();
  }

  void Block_Factor::factor(const Mpf_Batch &S, Mpf_Batch &Lf, Mpf_Batch &X,
                            Factor_Options options)
  {
    // n_ is a default only: one object per batch size can serve every block
    // size, so that its cached products are shared
    if(S.rows % 32 != 0 || S.cols != S.rows || S.batch != batch_)
      throw std::invalid_argument("Block_Factor::factor: shape");
    const int64_t n = S.rows;
    skip_triangular_zeros_ = options.skip_triangular_zeros;
    lower_update_ = options.lower_update;
    reuse_residues_ = options.reuse_residues;
    log_panels_ = options.log_panels;
    measure_residual_ = options.measure_residual;
    structured_solve_ = options.structured_solve;
    correction_limit_ = options.correction_limit;
    refine_min_ = options.refine_min_n;
    refine_max_ = options.refine_max_n;
    panels_skipped_ = 0;
    corrupt_leaf_limbs_ = options.corrupt_leaf_limbs;
    if(options.correction_limit > 0
       || (options.correction_limit < 0 && !options.equilibrate))
      throw std::invalid_argument(
        "correction_limit: negative log2, with global equilibration");
    if(options.correction_digits < 0 && !options.equilibrate)
      throw std::invalid_argument(
        "adaptive correction requires global equilibration");
    if(options.panel_refine < 0)
      throw std::invalid_argument("negative refinement count");
    if(options.symmetric_update_min < 0
       || (options.symmetric_update_min && batch_ != 1))
      throw std::invalid_argument("experimental symmetric update requires "
                                  "batch 1 and nonnegative threshold");
    if(options.correction_digits > D_
       || (options.correction_digits > 0 && options.correction_digits % 4))
      throw std::invalid_argument(
        "correction digits must be a multiple of four <= D");
    pool_active_ = true;
    product_guard_ = options.crt_guard;
    tic();
    mpf_zero(Lf, stream_);
    mpf_zero(X, stream_);
    toc("copy: zero");
    if(options.equilibrate)
      {
        Mpf_Batch A = S;
        int64_t *scale;
        BF_CHECK(cudaMallocAsync(&A.e, S.entries() * 8, stream_));
        BF_CHECK(cudaMallocAsync(&scale, int64_t(batch_) * n * 8, stream_));
        tic();
        k_equilibrate<<<nb(S.entries()), TH, 0, stream_>>>(S.e, A.e, scale, n,
                                                           batch_, L_);
        BF_CHECK(cudaGetLastError());
        toc("scale");
        const auto factor_scaled = [&] {
          rec(A, Lf, X, options.panel_refine, options.refine_min_n,
              options.correction_digits, options.symmetric_update_min);
        };
        if(corrects_always_.count(S.rows))
          {
            correction_limit_ = 0;
            refine_min_ = refine_max_ = 0;
          }
        try
          {
            factor_scaled();
          }
        catch(const Correction_Needed &)
          {
            // this matrix needs the correction: start over with every panel
            // corrected, now and from now on, for this shape
            ++correction_restarts_;
            corrects_always_.insert(S.rows);
            correction_limit_ = 0;
            refine_min_ = refine_max_ = 0;
            panels_skipped_ = 0;
            for(const auto &kv : live_)
              pool_put(kv.first, kv.second);
            live_.clear();
            factor_scaled();
          }
        tic();
        k_unscale<<<nb(S.entries()), TH, 0, stream_>>>(Lf.e, X.e, scale, n,
                                                       batch_);
        BF_CHECK(cudaGetLastError());
        toc("scale");
        BF_CHECK(cudaFreeAsync(A.e, stream_));
        BF_CHECK(cudaFreeAsync(scale, stream_));
      }
    else
      rec(S, Lf, X, options.panel_refine, options.refine_min_n,
          options.correction_digits, options.symmetric_update_min);
    pool_active_ = false;
    product_guard_ = 0;
    live_.clear();
    pool_flush();
    check();
  }

  // P = X B when X is lower triangular: columns of B with the same nonzero
  // rows (in units of 16 rows, over all blocks) form a group; the group's
  // product runs over those rows only, and the rows of P above the first one
  // are zero. Most of B is exact zeros in SDPB: the rows of a block with
  // several (r, s) pairs come in segments, and a column is nonzero on few
  // segments (zero entries of P are 36% at nmax14 and 18).
  void Block_Factor::structured_product(const Mpf_Batch &Xt, const Mpf_Batch &B,
                                        int32_t *d_m, int64_t *d_e, bool flat)
  {
    const int64_t buckets = (B.rows + 15) / 16;
    tic();
    uint8_t *d;
    std::vector<uint8_t> h(size_t(buckets) * B.cols);
    BF_CHECK(cudaMallocAsync(&d, h.size(), stream_));
    BF_CHECK(cudaMemsetAsync(d, 0, h.size(), stream_));
    k_column_rows<<<dim3(unsigned(B.cols), unsigned((B.rows + TH - 1) / TH),
                         unsigned(B.batch)),
                    TH, 0, stream_>>>(B.e, B.rows, B.cols, buckets, d);
    BF_CHECK(cudaGetLastError());
    BF_CHECK(cudaMemcpyAsync(h.data(), d, h.size(), cudaMemcpyDeviceToHost, stream_));
    BF_CHECK(cudaStreamSynchronize(stream_));
    BF_CHECK(cudaFreeAsync(d, stream_));
    // columns by the hull of their nonzero rows, in units of 16 rows
    struct Hull
    {
      int first, last; // rows first..last-1 of B, rows first.. of P
      std::vector<int> cols;
    };
    std::map<std::pair<int, int>, std::vector<int>> by_hull;
    for(int64_t j = 0; j < B.cols; ++j)
      {
        int first = -1, last = -1;
        for(int64_t q = 0; q < buckets; ++q)
          if(h[j * buckets + q])
            {
              if(first < 0)
                first = int(q);
              last = int(q) + 1;
            }
        if(first >= 0)
          by_hull[{first * 16, last * 16}].push_back(int(j));
      }
    std::vector<Hull> hulls;
    for(auto &kv : by_hull)
      hulls.push_back({kv.first.first, kv.first.second, std::move(kv.second)});
    // Merge hulls while that lowers the estimated time (in ms;
    // fit on real nmax14 blocks): a panel costs launches, the outputs
    // cost conversions and memory, the rows of B cost INT8 work and residues.
    const auto &basis = gemm_for(int(Xt.cols), int(B.cols), int(Xt.rows), D_).gemm();
    const double r = double(basis.basis().primes.size());
    const int np = basis.panel_width();
    const auto cost = [&](int first, int last, size_t cols) {
      const double panels = std::ceil(double(cols) / np), bn = double(batch_);
      const double groups = std::ceil(r / opt_.prime_group);
      const double launch = 0.02 + 1.2e-4 * 3 * opt_.prime_group * bn;
      const double rows = double(Xt.cols - first), inner = double(last - first);
      return panels * (groups * launch + 0.5)
             + bn * rows * cols * 2.7e-6
             + bn * rows * cols * inner * 6 * r / 2.3e12
             + bn * inner * cols * 0.5e-6;
    };
    while(hulls.size() > 1)
      {
        double best = 0;
        size_t a = 0, b = 0;
        for(size_t u = 0; u < hulls.size(); ++u)
          for(size_t v = u + 1; v < hulls.size(); ++v)
            {
              const double gain
                = cost(hulls[u].first, hulls[u].last, hulls[u].cols.size())
                  + cost(hulls[v].first, hulls[v].last, hulls[v].cols.size())
                  - cost(std::min(hulls[u].first, hulls[v].first),
                         std::max(hulls[u].last, hulls[v].last),
                         hulls[u].cols.size() + hulls[v].cols.size());
              if(gain > best)
                {
                  best = gain;
                  a = u;
                  b = v;
                }
            }
        if(best <= 0)
          break;
        hulls[a].first = std::min(hulls[a].first, hulls[b].first);
        hulls[a].last = std::max(hulls[a].last, hulls[b].last);
        hulls[a].cols.insert(hulls[a].cols.end(), hulls[b].cols.begin(),
                             hulls[b].cols.end());
        hulls.erase(hulls.begin() + b);
      }
    std::vector<Mpf_Gemm::Column_Group> groups;
    for(auto &u : hulls)
      {
        Mpf_Gemm::Column_Group g;
        g.cols = std::move(u.cols);
        g.first = u.first;
        g.last = u.last;
        groups.push_back(std::move(g));
      }
    toc("solve: column groups");
    if(opt_.timing)
      ms_["count/column groups"] += groups.size();
    // the zeros of P (written into the planes, whichever the layout)
    tic();
    {
      Mpf_Batch Z;
      Z.m = d_m;
      Z.e = d_e;
      Z.L = L_;
      Z.batch = batch_;
      Z.rows = Xt.cols;
      Z.cols = B.cols;
      mpf_zero(Z, stream_);
    }
    toc("solve: zero");
    tic();
    auto &g = gemm_for(int(Xt.cols), int(B.cols), int(Xt.rows), D_);
    std::map<std::string, double> before, before2;
    if(opt_.profile_products)
      {
        before = g.timings();
        before2 = g.gemm().timings();
      }
    g.multiply_groups(Xt.view(), B.view(), groups, d_m, d_e, flat);
    toc("role/" + role_ + "/k=" + std::to_string(Xt.rows) + " groups");
    if(opt_.timing)
      ms_["count/" + role_ + "/groups"] += 1;
    if(opt_.profile_products)
      {
        for(const auto &kv : g.timings())
          ms_["kernel/" + role_ + "/" + kv.first] += kv.second - before[kv.first];
        for(const auto &kv : g.gemm().timings())
          ms_["kernel/" + role_ + "/" + kv.first] += kv.second - before2[kv.first];
      }
  }

  void Block_Factor::solve(const Mpf_Batch &Lf, const Mpf_Batch &X,
                           const Mpf_Batch &B, Mpf_Batch &P, int refine)
  {
    role_ = "solve P=X B";
    if(structured_solve_)
      {
        Mpf_Batch Xt = transposed(X);
        if(!P.m)
          P = alloc(X.rows, B.cols);
        structured_product(Xt, B, P.m, P.e, false);
        release(Xt);
      }
    else
      mm(X, false, B, false, P, 0, 1);
    for(int it = 0; it < refine; ++it)
      {
        Mpf_Batch LP, Rm = alloc(B.rows, B.cols), D;
        mm(Lf, false, P, false, LP);
        mpf_batch_add(B, LP, -1, Rm, stream_); // R = B - L P
        release(LP);
        mm(X, false, Rm, false, D);
        mpf_batch_add(P, D, 1, P, stream_); // P += X R
        release(Rm);
        release(D);
      }
    check();
  }
}
