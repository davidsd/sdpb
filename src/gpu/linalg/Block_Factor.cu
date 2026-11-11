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
#include "../rns/rns_kernels.cuh"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

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

    __global__ void
    k_transpose(const int32_t *__restrict__ am, const int64_t *__restrict__ ae,
                int L, int batch, int64_t R, int64_t C,
                int32_t *__restrict__ tm, int64_t *__restrict__ te)
    {
      const int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      const int64_t tot = int64_t(batch) * R * C;
      if(idx >= tot)
        return;
      const int64_t b = idx / (R * C), rem = idx % (R * C);
      const int64_t i = rem % R, j = rem / R;    // A(i, j)
      const int64_t out = b * R * C + i * C + j; // T(j, i), T has C rows
      te[out] = ae[idx];
      for(int t = 0; t < L; ++t)
        tm[(int64_t)t * tot + out] = am[(int64_t)t * tot + idx];
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
    __global__ void
    k_add(const int32_t *__restrict__ Am, const int64_t *__restrict__ Ae,
          const int32_t *__restrict__ Bm, const int64_t *__restrict__ Be,
          int sgn, int L, int64_t n, int32_t *__restrict__ Om,
          int64_t *__restrict__ Oe)
    {
      const int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(i >= n)
        return;
      uint32_t a[LM], b[LM], wa[LM + 2], wb[LM + 2];
      for(int k = 0; k < L; ++k)
        {
          a[k] = uint32_t(Am[(int64_t)k * n + i]);
          b[k] = uint32_t(Bm[(int64_t)k * n + i]);
        }
      bool na = a[L - 1] >> 31, nbg = ((b[L - 1] >> 31) != 0) ^ (sgn < 0);
      if(a[L - 1] >> 31)
        negate(a, L);
      if(b[L - 1] >> 31)
        negate(b, L);
      long long ea = Ae[i], eb = Be[i];
      if(bitlen(a, L) == 0)
        ea = MPF_EZ;
      if(bitlen(b, L) == 0)
        eb = MPF_EZ;
      if(ea == MPF_EZ && eb == MPF_EZ)
        {
          for(int k = 0; k < L; ++k)
            Om[(int64_t)k * n + i] = 0;
          Oe[i] = MPF_EZ;
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
      store_norm(wa, L + 2, t, L, Om, n, i, Oe);
    }

    __global__ void k_zero_e(int64_t *e, int64_t n)
    {
      const int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(i < n)
        e[i] = MPF_EZ;
    }

    __global__ void
    k_max_exp(const int64_t *e, int64_t count, unsigned long long *out)
    {
      __shared__ long long v[256];
      long long mx = LLONG_MIN;
      for(int64_t q = int64_t(blockIdx.x) * blockDim.x + threadIdx.x;
          q < count; q += int64_t(gridDim.x) * blockDim.x)
        if(e[q] != MPF_EZ)
          mx = max(mx, (long long)e[q]);
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
              int L, int tiles, int WL, int LF, uint32_t *__restrict__ W,
              int64_t *__restrict__ kk)
    {
      const int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      const int64_t tot = int64_t(tiles) * 1024;
      if(idx >= tot)
        return;
      const int64_t b = idx / 1024, rem = idx % 1024;
      const int i = int(rem % 32), j = int(rem / 32); // column-major input
      auto kdiag = [&](int r) -> long long {
        const int64_t d = b * 1024 + r * 32 + r;
        // normalized: |m| has 32L-2 bits; a_rr < 2^(e + 32L - 2)
        long long top = e[d] + 32LL * L - 2;
        return (top + 1) >> 1; // ceil(top / 2), floor for negatives
      };
      const long long ki = kdiag(i), kj = kdiag(j);
      if(j == 0)
        kk[b * 32 + i] = ki;
      uint32_t a[LM], w[WM];
      for(int t = 0; t < L; ++t)
        a[t] = uint32_t(m[(int64_t)t * tot + idx]);
      uint32_t *out = W + ((b * 32 + i) * 32 + j) * WL;
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
               int32_t *__restrict__ m, int64_t *__restrict__ e)
    {
      const int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      const int64_t tot = int64_t(tiles) * 1024;
      if(idx >= tot)
        return;
      const int64_t b = idx / 1024, rem = idx % 1024;
      const int i = int(rem % 32), j = int(rem / 32);
      uint32_t w[WM];
      if(i < j)
        {
          for(int t = 0; t < L; ++t)
            m[(int64_t)t * tot + idx] = 0;
          e[idx] = MPF_EZ;
          return;
        }
      const uint32_t *in = W + ((b * 32 + i) * 32 + j) * WL;
      for(int t = 0; t < WL; ++t)
        w[t] = in[t];
      const long long sc = inverse ? -kk[b * 32 + j] : kk[b * 32 + i];
      store_norm(w, WL, -32LL * LF + sc, L, m, tot, idx, e);
    }
  }

  // ---------------------------------------------------------------- helpers
  void mpf_transpose(const Mpf_Batch &A, Mpf_Batch &T, cudaStream_t s)
  {
    k_transpose<<<nb(A.entries()), TH, 0, s>>>(A.m, A.e, A.L, A.batch, A.rows,
                                               A.cols, T.m, T.e);
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
                     Mpf_Batch &C, cudaStream_t s)
  {
    k_add<<<nb(A.entries()), TH, 0, s>>>(A.m, A.e, B.m, B.e, sgn, A.L,
                                         A.entries(), C.m, C.e);
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
    if(L != 28)
      throw std::invalid_argument(
        "Block_Factor: leaves instantiated for L = 28");
    BF_CHECK(cudaEventCreate(&ev0_));
    BF_CHECK(cudaEventCreate(&ev1_));
  }

  Block_Factor::~Block_Factor()
  {
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

  Mpf_Batch Block_Factor::alloc(int64_t rows, int64_t cols)
  {
    Mpf_Batch M;
    M.L = L_;
    M.batch = batch_;
    M.rows = rows;
    M.cols = cols;
    auto allocate = [&](void **p, size_t bytes) {
      auto error = malloc_async_retry(p, bytes, stream_);
      if(error == cudaErrorMemoryAllocation)
        {
          // No cached plan is borrowed across alloc(). Wait for all queued
          // products before releasing their buffers, then retry the temporary
          // allocation. Product creation already has a separate memory budget.
          (void)cudaGetLastError();
          check();
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
    BF_CHECK(allocate(reinterpret_cast<void **>(&M.m),
                       size_t(L_) * M.entries() * 4));
    auto error = allocate(reinterpret_cast<void **>(&M.e),
                            size_t(M.entries()) * 8);
    if(error != cudaSuccess)
      {
        BF_CHECK(cudaFreeAsync(M.m, stream_));
        BF_CHECK(error);
      }
    return M;
  }
  void Block_Factor::release(Mpf_Batch &M)
  {
    if(M.m)
      BF_CHECK(cudaFreeAsync(M.m, stream_));
    if(M.e)
      BF_CHECK(cudaFreeAsync(M.e, stream_));
    M.m = nullptr;
    M.e = nullptr;
  }

  Mpf_Gemm &Block_Factor::gemm_for(int m, int n, int k, int digits)
  {
    auto key = std::make_tuple(m, n, k, digits);
    auto it = gemms_.find(key);
    if(it == gemms_.end())
      {
        Gram_Options o = opt_;
        o.timing = false;
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
        BF_CHECK(cudaStreamSynchronize(stream_));
        // A caller may retain an async allocation pool. Its unused reserve is
        // not available to cudaMalloc and otherwise looks like live memory to
        // the planner, causing avoidable eviction of reusable product plans.
        int device;
        cudaMemPool_t pool;
        BF_CHECK(cudaGetDevice(&device));
        BF_CHECK(cudaDeviceGetDefaultMemPool(&pool, device));
        BF_CHECK(cudaMemPoolTrimTo(pool, 0));
        const double r = 132, g = o.prime_group, bm = double(batch_) * m;
        // residues of C: uint16 per prime, or with reduce_to_crt (groups of
        // at most 32 primes) the CRT rows (2 x 144 bytes) and a float
        const bool fused
          = o.reduce_to_crt && !o.fused_crt && o.prime_group <= 32;
        auto need = [&](int pn, bool ca) {
          return 2.0 * g * o.max_entries_per_residue_batch * 4 // residue batch
                 + 3.0 * g * bm * k
                 + 3.0 * g * batch_ * double(k) * pn   // RA, RB
                 + 12.0 * g * bm * pn                  // int32 products
                 + (fused ? 292.0 : 2.0 * r) * bm * pn // residues of C
                 + double(o.crt_batch) * (2 * 144 + 2 * 4 * 264) // CRT buffers
                 + 4.0 * 60 * bm * pn                            // exact panel
                 + (ca ? 3.0 * r * bm * k : 0.0); // cached residues of A
        };
        bool cache_a = o.cache_a;
        int panel = std::min(n, 1024);
        for(int pass = 0; pass < 2; ++pass)
          {
            size_t free_b = 0, total_b = 0;
            BF_CHECK(cudaMemGetInfo(&free_b, &total_b));
            // fraction of the free memory one product may plan for
            // (SDPB_GPU_BUDGET, default 0.7: on the real nmax14 blocks 0.6
            // left the residue cache off and cost 1.6x in P = X B; with all 36
            // nmax18 blocks at once 0.8 ran out of memory)
            const char *bf = getenv("SDPB_GPU_BUDGET");
            const double budget = (bf ? atof(bf) : 0.7) * double(free_b);
            cache_a = o.cache_a;
            panel = std::min(n, 1024);
            while(need(panel, cache_a) > budget && panel > 128)
              panel /= 2;
            if(need(panel, cache_a) > budget && cache_a)
              {
                cache_a = false;
                panel = std::min(n, 1024);
                while(need(panel, cache_a) > budget && panel > 128)
                  panel /= 2;
              }
            if(getenv("SDPB_GPU_VERBOSE"))
              printf("  [Block_Factor] product m=%d n=%d k=%d batch=%d: panel "
                     "%d, cache_a %d, "
                     "buffers %.2f GB, free %.2f GB, %zu cached products\n",
                     m, n, k, batch_, panel, int(cache_a),
                     need(panel, cache_a) / 1e9, free_b / 1e9, gemms_.size());
            if(need(panel, cache_a) <= budget || gemms_.empty())
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
    const int64_t m = tU ? U.cols : U.rows, k = tU ? U.rows : U.cols;
    const int64_t n = tV ? V.rows : V.cols;
    Mpf_Batch Ap = U, Bp = V;
    bool freeA = false, freeB = false;
    if(!tU)
      {
        Ap = alloc(U.cols, U.rows);
        mpf_transpose(U, Ap, stream_);
        freeA = true;
      }
    if(tV)
      {
        Bp = alloc(V.cols, V.rows);
        mpf_transpose(V, Bp, stream_);
        freeB = true;
      }
    if(!C.m)
      C = alloc(m, n);
    tic();
    gemm_for(int(m), int(n), int(k), digits)
      .multiply(Ap.view(), Bp.view(), C.m, C.e,
                skip_triangular_zeros_ ? triangles : 0);
    toc(std::string("product, m=") + std::to_string(m)
        + " n=" + std::to_string(n) + " k=" + std::to_string(k));
    if(freeA)
      release(Ap);
    if(freeB)
      release(Bp);
  }

  void Block_Factor::leaf(const Mpf_Batch &A, Mpf_Batch &Lf, Mpf_Batch &X)
  {
    constexpr int LF = 28, WL = LF + 24;
    const int tiles = A.batch;
    uint32_t *W, *R, *Xw, *S;
    int64_t *kk;
    const size_t wb = size_t(tiles) * 1024 * WL * 4;
    BF_CHECK(cudaMallocAsync(&W, wb, stream_));
    BF_CHECK(cudaMallocAsync(&Xw, wb, stream_));
    BF_CHECK(cudaMallocAsync(&S, wb, stream_));
    BF_CHECK(cudaMallocAsync(&R, size_t(tiles) * 32 * WL * 4, stream_));
    BF_CHECK(cudaMallocAsync(&kk, size_t(tiles) * 32 * 8, stream_));
    tic();
    k_leaf_in<<<nb(int64_t(tiles) * 1024), TH, 0, stream_>>>(
      A.m, A.e, L_, tiles, WL, LF, W, kk);
    BF_CHECK(cudaGetLastError());
    if(opt_.fused_leaf && opt_.leaf_lookahead)
      leaf_chol_inv<WL, LF, 32, 512, true, true>
        <<<tiles, 512, 0, stream_>>>(W, Xw, nullptr);
    else if(opt_.fused_leaf)
      leaf_chol_inv<WL, LF, 32, 512, true, false>
        <<<tiles, 512, 0, stream_>>>(W, Xw, nullptr);
    else
      {
        chol_tiles<WL, LF, 32><<<tiles, 256, 0, stream_>>>(W, R);
        BF_CHECK(cudaGetLastError());
        trinv_tiles<WL, LF, 32><<<tiles, 256, 0, stream_>>>(W, R, Xw, S);
      }
    BF_CHECK(cudaGetLastError());
    k_leaf_out<<<nb(int64_t(tiles) * 1024), TH, 0, stream_>>>(
      W, WL, LF, kk, tiles, false, L_, Lf.m, Lf.e);
    k_leaf_out<<<nb(int64_t(tiles) * 1024), TH, 0, stream_>>>(
      Xw, WL, LF, kk, tiles, true, L_, X.m, X.e);
    BF_CHECK(cudaGetLastError());
    for(void *p : {(void *)W, (void *)Xw, (void *)S, (void *)R, (void *)kk})
      BF_CHECK(cudaFreeAsync(p, stream_));
    toc("leaves (32 x 32 Cholesky + inverse)");
  }

  void Block_Factor::refined_panel(const Mpf_Batch &B, const Mpf_Batch &Lf,
                                   const Mpf_Batch &X, Mpf_Batch &Z,
                                   int refine, int cdigits)
  {
    mm(B, false, X, true, Z, 0, 4);
    for(int it = 0; it < refine; ++it)
      {
        Mpf_Batch ZL, R = alloc(B.rows, B.cols), C;
        mm(Z, false, Lf, true, ZL, 0, 4);
        mpf_batch_add(B, ZL, -1, R, stream_);
        mm(R, false, X, true, C,
           cdigits < 0 ? correction_digits(R, X) : cdigits, 4);
        mpf_batch_add(Z, C, 1, Z, stream_);
        release(ZL);
        release(R);
        release(C);
      }
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

  // (Lf, X) of A (n x n): SDPviBe's chol_inv recursion
  void Block_Factor::rec(const Mpf_Batch &A, Mpf_Batch &Lf, Mpf_Batch &X,
                         int refine, int min_n, int cdigits, int gram_min)
  {
    const int64_t n = A.rows;
    if(n == 32)
      {
        leaf(A, Lf, X);
        return;
      }
    const int64_t n1 = (n / 2 + 31) / 32 * 32, n2 = n - n1;
    Mpf_Batch A11 = alloc(n1, n1), A21 = alloc(n2, n1), A22 = alloc(n2, n2);
    mpf_sub_block(A, 0, 0, A11, false, stream_);
    mpf_sub_block(A, n1, 0, A21, false, stream_);
    mpf_sub_block(A, n1, n1, A22, false, stream_);
    Mpf_Batch L11 = alloc(n1, n1), X11 = alloc(n1, n1);
    rec(A11, L11, X11, refine, min_n, cdigits, gram_min);
    release(A11);
    Mpf_Batch L21, Wm, T, X21;
    refined_panel(A21, L11, X11, L21, n1 >= min_n ? refine : 0, cdigits);
    release(A21);
    if(gram_min && n2 >= gram_min)
      symmetric_product(L21, Wm);
    else
      mm(L21, false, L21, true, Wm); // L21 L21^T
    Mpf_Batch Sm = alloc(n2, n2);
    mpf_batch_add(A22, Wm, -1, Sm, stream_); // S = A22 - L21 L21^T
    release(A22);
    release(Wm);
    Mpf_Batch L22 = alloc(n2, n2), X22 = alloc(n2, n2);
    rec(Sm, L22, X22, refine, min_n, cdigits, gram_min);
    release(Sm);
    mm(L21, false, X11, false, T, 0, 2); // L21 X11
    mm(X22, false, T, false, X21, 0, 1); // X22 L21 X11
    release(T);
    mpf_zero(Lf, stream_);
    mpf_zero(X, stream_);
    mpf_sub_block(Lf, 0, 0, L11, true, stream_);
    mpf_sub_block(Lf, n1, 0, L21, true, stream_);
    mpf_sub_block(Lf, n1, n1, L22, true, stream_);
    Mpf_Batch Z = alloc(n2, n1);
    mpf_zero(Z, stream_);
    mpf_batch_add(Z, X21, -1, X21, stream_); // X21 = -X22 L21 X11
    release(Z);
    mpf_sub_block(X, 0, 0, X11, true, stream_);
    mpf_sub_block(X, n1, 0, X21, true, stream_);
    mpf_sub_block(X, n1, n1, X22, true, stream_);
    for(Mpf_Batch *M : {&L11, &X11, &L21, &L22, &X22, &X21})
      release(*M);
  }

  // Congruence-scaled inputs have diagonal in [1/4,1), hence every exact
  // Cholesky row has norm < 1. Budget the absolute error of the SMALL correction
  // below 2^-900: k*64*2^(-8D)*max|R|*max|X| <= 2^-900.
  // The full-precision residual is untouched. Cap D at the available precision;
  // reaching that cap does not certify arbitrary ill-conditioned inputs.
  int Block_Factor::correction_digits(const Mpf_Batch &R, const Mpf_Batch &X)
  {
    unsigned long long *d, h[2] = {0, 0};
    BF_CHECK(cudaMallocAsync(&d, 16, stream_));
    BF_CHECK(cudaMemsetAsync(d, 0, 16, stream_));
    k_max_exp<<<256, 256, 0, stream_>>>(R.e, R.entries(), d);
    k_max_exp<<<256, 256, 0, stream_>>>(X.e, X.entries(), d + 1);
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
                       : 900 + 6 + bits_k + er + ex + 2 * (32 * L_ - 2);
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
    if(S.rows != n_ || S.cols != n_ || S.batch != batch_)
      throw std::invalid_argument("Block_Factor::factor: shape");
    skip_triangular_zeros_ = options.skip_triangular_zeros;
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
    if(options.equilibrate)
      {
        Mpf_Batch A = S;
        int64_t *scale;
        BF_CHECK(cudaMallocAsync(&A.e, S.entries() * 8, stream_));
        BF_CHECK(cudaMallocAsync(&scale, int64_t(batch_) * n_ * 8, stream_));
        k_equilibrate<<<nb(S.entries()), TH, 0, stream_>>>(S.e, A.e, scale, n_,
                                                           batch_, L_);
        BF_CHECK(cudaGetLastError());
        rec(A, Lf, X, options.panel_refine, options.refine_min_n,
            options.correction_digits, options.symmetric_update_min);
        k_unscale<<<nb(S.entries()), TH, 0, stream_>>>(Lf.e, X.e, scale, n_,
                                                       batch_);
        BF_CHECK(cudaGetLastError());
        BF_CHECK(cudaFreeAsync(A.e, stream_));
        BF_CHECK(cudaFreeAsync(scale, stream_));
      }
    else
      rec(S, Lf, X, options.panel_refine, options.refine_min_n,
          options.correction_digits, options.symmetric_update_min);
    check();
  }

  void Block_Factor::solve(const Mpf_Batch &Lf, const Mpf_Batch &X,
                           const Mpf_Batch &B, Mpf_Batch &P, int refine)
  {
    mm(X, false, B, false, P);
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
