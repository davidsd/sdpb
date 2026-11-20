// GPU backend for SDPB: Q-formation path on device floats.  See Mpf_Gram.hxx.
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: the float format, the column alignment and the
// rounding back to floats follow sdpvibe/mpf.py (MPF, normalize, GramAcc) and
// sdpvibe/limbops.py (shift_round).  Portions adapted from SDPviBe, Copyright
// (c) 2026 Antoine Tilloy, MIT License (docs/GPU_backend_credits.md).

#include "Mpf_Gram.hxx"
#include "rns_kernels.cuh"

#include <climits>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdexcept>

#define MG_CHECK(x)                                                           \
  do                                                                          \
    {                                                                         \
      cudaError_t e_ = (x);                                                   \
      if(e_ != cudaSuccess)                                                   \
        {                                                                     \
          std::ostringstream ss_;                                             \
          ss_ << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(e_); \
          throw std::runtime_error(ss_.str());                                \
        }                                                                     \
  } while(0)

namespace sdpb_gpu
{
  namespace
  {
    constexpr int threads = 256;
    constexpr int LMAX = 64; // limbs per float handled by the kernels

    int64_t nblocks(int64_t n) { return (n + threads - 1) / threads; }

    // bit length of |x| for an L-limb two's complement x (|x| = 2^k for
    // negatives may read one less; callers add a guard bit)
    __device__ int bitlen_tc(const int32_t *m, int64_t plane, int64_t idx,
                             int L)
    {
      uint32_t ext = m[(int64_t)(L - 1) * plane + idx] < 0 ? 0xffffffffu : 0u;
      for(int t = L - 1; t >= 0; --t)
        {
          uint32_t w = uint32_t(m[(int64_t)t * plane + idx]) ^ ext;
          if(w)
            return 32 * t + 32 - __clz(w);
        }
      return 0;
    }

    // column maxima of e + bitlen(m) (+1 guard bit): one thread block per
    // (column, slice of rows), a block reduction, then one atomicMax
    __global__ void col_tops(const int32_t *__restrict__ m,
                             const int64_t *__restrict__ e, int L,
                             int64_t rows, int64_t cols, int64_t pl,
                             unsigned long long *top, int64_t in_bstride = 0,
                             int64_t top_bstride = 0)
    {
      // batch b = blockIdx.z: block b's matrix at offset b*in_bstride, its
      // column maxima at top + b*top_bstride
      m += blockIdx.z * in_bstride;
      e += blockIdx.z * in_bstride;
      top += blockIdx.z * top_bstride;
      const int64_t j = blockIdx.x;
      const int64_t r0 = int64_t(blockIdx.y) * blockDim.x * 16;
      long long best = LLONG_MIN;
      for(int64_t i = r0 + threadIdx.x; i < rows && i < r0 + blockDim.x * 16;
          i += blockDim.x)
        {
          const int64_t idx = j * rows + i;
          const int64_t ex = e[idx];
          if(ex == MPF_EZ)
            continue;
          int bl = bitlen_tc(m, pl, idx, L);
          if(bl)
            best = max(best, (long long)(ex + bl + 1));
        }
      __shared__ long long red[256];
      red[threadIdx.x] = best;
      __syncthreads();
      for(int w = blockDim.x / 2; w > 0; w /= 2)
        {
          if(threadIdx.x < w)
            red[threadIdx.x] = max(red[threadIdx.x], red[threadIdx.x + w]);
          __syncthreads();
        }
      if(threadIdx.x == 0 && red[0] != LLONG_MIN)
        atomicMax(top + j, (unsigned long long)(red[0]) ^ 0x8000000000000000ULL);
    }

    __global__ void make_refs(const unsigned long long *top, int64_t *Eref,
                              int64_t cols, int D)
    {
      int64_t j = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(j >= cols)
        return;
      long long t = (long long)(top[j] ^ 0x8000000000000000ULL);
      Eref[j] = t <= MPF_EZ / 2 ? 0 : t - 8LL * D + 1;
    }

    // Residues of A = round_half_up(m 2^(e - Eref_j)) (|A| < 2^(8D-2)) modulo
    // gc primes, straight from the float's limbs, all in registers:
    //   w = m sign-extended to L+1 limbs (+ 2^(s-1) if s >= 1), shifted right by
    //   s = 32 qs + rs bits with a barrel shifter over limbs (static indices,
    //   predicated) and funnel shifts, giving A as WA = L+2 two's complement
    //   limbs; then res(A) = sum_k a_k c_k - neg c_WA (mod m), c_k = 2^(32k)
    //   mod m from shared memory with the same index in every lane (broadcast).
    // Centred, r = 128 h + l, s = h + l, written to R[(X*gc + p)*plane +
    // j*ld + i]; rows i >= rows are zero.
    template <int LT>
    __global__ void mpf_residues(const int32_t *__restrict__ m,
                                 const int64_t *__restrict__ e, int Lrt,
                                 int64_t rows, int64_t cols, int64_t ld,
                                 int64_t pl, int64_t r_off,
                                 const int64_t *__restrict__ Eref,
                                 const uint32_t *__restrict__ ctab,
                                 const uint32_t *__restrict__ unused,
                                 const int32_t *__restrict__ primes,
                                 const double *__restrict__ minv, int gc,
                                 int8_t *__restrict__ R, int64_t plane,
                                 int64_t in_bstride = 0, int64_t eref_bstride = 0,
                                 int64_t r_bstride = 0)
    {
      // batch b = blockIdx.y: input at offset b*in_bstride, column exponents
      // at Eref + b*eref_bstride, output at r_off + b*r_bstride
      m += blockIdx.y * in_bstride;
      e += blockIdx.y * in_bstride;
      Eref += blockIdx.y * eref_bstride;
      r_off += blockIdx.y * r_bstride;
      const int L = LT > 0 ? LT : Lrt;
      constexpr int XN = LT > 0 ? LT + 3 : LMAX + 3; // working limbs
      const int CK = L + 3, WA = L + 2;
      extern __shared__ uint32_t sc[]; // [gc][CK]
      for(int q = threadIdx.x; q < gc * CK; q += blockDim.x)
        sc[q] = ctab[q];
      __syncthreads();
      const int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(idx >= ld * cols)
        return;
      const int64_t j = idx / ld, i = idx % ld;
      uint32_t x[XN];
      bool zero = true;
      uint32_t ext = 0;
      if(i < rows)
        {
          const int64_t in = j * rows + i;
          const int64_t ex = e[in];
          const long long s = (ex == MPF_EZ) ? (1LL << 40) : (long long)(Eref[j] - ex);
          if(s <= 32LL * L + 1 && s >= -32LL)
            {
              zero = false;
              // x = [0, w_0, ..., w_L] (one limb up, so that s >= -32 is a right shift)
#pragma unroll
              for(int t = 0; t < XN; ++t)
                x[t] = (t >= 1 && t <= L) ? uint32_t(m[(int64_t)(t - 1) * pl + in]) : 0u;
              ext = int32_t(x[L]) < 0 ? 0xffffffffu : 0u;
#pragma unroll
              for(int t = 0; t < XN; ++t)
                if(t >= L + 1)
                  x[t] = ext;
              // rounding: + 2^(s-1) at bit s - 1 + 32 of x
              if(s >= 1)
                {
                  const int qr = int((s + 31) >> 5), rr = int((s + 31) & 31);
                  unsigned long long c = 0;
#pragma unroll
                  for(int t = 0; t < XN; ++t)
                    {
                      unsigned long long add = (t == qr) ? (1ULL << rr) : 0ULL;
                      unsigned long long v = (unsigned long long)x[t] + add + c;
                      x[t] = uint32_t(v);
                      c = v >> 32;
                    }
                }
              ext = int32_t(x[XN - 1]) < 0 ? 0xffffffffu : 0u;
              // right shift by s + 32 = 32 u + rs bits, u in [0, L+1]
              const int u = int((s + 32) >> 5), rs = int((s + 32) & 31);
#pragma unroll
              for(int b = 0; b < 7; ++b)
                {
                  const int st = 1 << b;
                  if(st < XN && (u & st))
                    {
#pragma unroll
                      for(int t = 0; t < XN; ++t)
                        x[t] = (t + st < XN) ? x[t + st] : ext;
                    }
                }
#pragma unroll
              for(int t = 0; t < XN; ++t)
                x[t] = __funnelshift_r(x[t], (t + 1 < XN) ? x[t + 1] : ext, rs);
            }
        }
      const bool neg = int32_t(x[WA - 1]) < 0 && !zero;
      for(int p = 0; p < gc; ++p)
        {
          int8_t hi8 = 0, lo8 = 0;
          if(!zero)
            {
              const int32_t mp = primes[p];
              const uint32_t *c = sc + p * CK;
              unsigned long long acc = neg ? (unsigned long long)(mp - c[WA]) : 0ULL;
#pragma unroll
              for(int t = 0; t < XN; ++t)
                if(t < WA)
                  acc += (unsigned long long)x[t] * c[t];
              int32_t r = rns_detail::mod_d((long long)acc, mp, minv[p]);
              if(r > mp / 2)
                r -= mp;
              int32_t lo = ((r + 64) & 127) - 64;
              hi8 = (int8_t)((r - lo) >> 7);
              lo8 = (int8_t)lo;
            }
          R[(int64_t)(0 * gc + p) * plane + r_off + idx] = hi8;
          R[(int64_t)(1 * gc + p) * plane + r_off + idx] = lo8;
          R[(int64_t)(2 * gc + p) * plane + r_off + idx] = (int8_t)(hi8 + lo8);
        }
    }

    // exact packed Q (W limb planes) * 2^(Eref_i + Eref_j) -> floats with
    // |m| of exactly 32L-2 bits (round half away from zero)
    __global__ void normalize_packed(const int32_t *__restrict__ Qi, int W,
                                     int64_t total, int64_t N,
                                     const int64_t *__restrict__ Eref, int L,
                                     int32_t *__restrict__ om,
                                     int64_t *__restrict__ oe)
    {
      int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(e >= total)
        return;
      // column j of packed entry e: j N - j(j-1)/2 <= e
      double b = 2.0 * N + 1.0;
      int64_t j = (int64_t)floor((b - sqrt(b * b - 8.0 * double(e))) / 2.0);
      auto off = [&](int64_t c) { return c * N - c * (c - 1) / 2; };
      while(j > 0 && off(j) > e)
        --j;
      while(off(j + 1) <= e)
        ++j;
      const int64_t i = j + (e - off(j));
      constexpr int WMAX = 96;
      uint32_t a[WMAX];
      for(int t = 0; t < W; ++t)
        a[t] = uint32_t(Qi[(int64_t)t * total + e]);
      const bool neg = int32_t(a[W - 1]) < 0;
      if(neg)
        {
          unsigned long long c = 1;
          for(int t = 0; t < W; ++t)
            {
              unsigned long long v = (unsigned long long)(~a[t]) + c;
              a[t] = uint32_t(v);
              c = v >> 32;
            }
        }
      int bl = 0;
      for(int t = W - 1; t >= 0; --t)
        if(a[t])
          {
            bl = 32 * t + 32 - __clz(a[t]);
            break;
          }
      if(bl == 0)
        {
          for(int t = 0; t < L; ++t)
            om[(int64_t)t * total + e] = 0;
          oe[e] = MPF_EZ;
          return;
        }
      const int NB = 32 * L - 2;
      int sh = bl - NB;
      uint32_t r[LMAX + 1];
      for(int pass = 0; pass < 2; ++pass)
        {
          // r = round_half_up(a / 2^sh) on the magnitude (sh may be < 0)
          for(int t = 0; t <= L; ++t)
            {
              long long bit = 32LL * t + sh;
              long long q = bit >> 5;
              int rr = int(bit - (q << 5));
              uint32_t lo = (q >= 0 && q < W) ? a[q] : 0u;
              uint32_t hi = (q + 1 >= 0 && q + 1 < W) ? a[q + 1] : 0u;
              r[t] = uint32_t((((unsigned long long)hi << 32) | lo) >> rr);
            }
          if(sh >= 1)
            {
              long long bit = sh - 1;
              uint32_t g = (bit >> 5) < W ? (a[bit >> 5] >> (bit & 31)) & 1u : 0u;
              unsigned long long c = g;
              for(int t = 0; t <= L && c; ++t)
                {
                  unsigned long long v = (unsigned long long)r[t] + c;
                  r[t] = uint32_t(v);
                  c = v >> 32;
                }
            }
          // renormalize once if rounding carried into bit NB
          if((r[L - 1] >> (NB - 32 * (L - 1))) == 0 && r[L] == 0)
            break;
          ++sh;
        }
      if(neg)
        {
          unsigned long long c = 1;
          for(int t = 0; t < L; ++t)
            {
              unsigned long long v = (unsigned long long)(~r[t]) + c;
              r[t] = uint32_t(v);
              c = v >> 32;
            }
        }
      for(int t = 0; t < L; ++t)
        om[(int64_t)t * total + e] = int32_t(r[t]);
      oe[e] = Eref[i] + Eref[j] + sh;
    }
    // exact block products (W limb planes over batch*mp*n, column-major per
    // block, leading dimension mp) times 2^(EA_i + EB_j) -> floats (batch, m x n)
    // panel version: columns j0..j0+pc-1, exact input (b, jj, i) at
    // b*mp*pc + jj*mp + i of in_plane entries; output (batch, m x n) floats
    // one entry of normalize_blocks: exact limbs Ci[t*in_plane + ei] (W) times
    // 2^ebase -> float om[t*total + eo], oe[eo]
    __device__ __noinline__ void round_entry(const int32_t *__restrict__ Ci,
                                             int W, int64_t in_plane,
                                             int64_t ei, int L, int64_t total,
                                             int64_t eo, int64_t ebase,
                                             int32_t *__restrict__ om,
                                             int64_t *__restrict__ oe)
    {
      constexpr int WMAX = 96;
      uint32_t a[WMAX];
      for(int t = 0; t < W; ++t)
        a[t] = uint32_t(Ci[(int64_t)t * in_plane + ei]);
      const bool neg = int32_t(a[W - 1]) < 0;
      if(neg)
        {
          unsigned long long c = 1;
          for(int t = 0; t < W; ++t)
            {
              unsigned long long v = (unsigned long long)(~a[t]) + c;
              a[t] = uint32_t(v);
              c = v >> 32;
            }
        }
      int bl = 0;
      for(int t = W - 1; t >= 0; --t)
        if(a[t])
          {
            bl = 32 * t + 32 - __clz(a[t]);
            break;
          }
      if(bl == 0)
        {
          for(int t = 0; t < L; ++t)
            om[(int64_t)t * total + eo] = 0;
          oe[eo] = MPF_EZ;
          return;
        }
      const int NB = 32 * L - 2;
      int sh = bl - NB;
      uint32_t r[LMAX + 1];
      for(int pass = 0; pass < 2; ++pass)
        {
          for(int t = 0; t <= L; ++t)
            {
              long long bit = 32LL * t + sh;
              long long q = bit >> 5;
              int rr = int(bit - (q << 5));
              uint32_t lo = (q >= 0 && q < W) ? a[q] : 0u;
              uint32_t hi = (q + 1 >= 0 && q + 1 < W) ? a[q + 1] : 0u;
              r[t] = uint32_t((((unsigned long long)hi << 32) | lo) >> rr);
            }
          if(sh >= 1)
            {
              long long bit = sh - 1;
              uint32_t g = (bit >> 5) < W ? (a[bit >> 5] >> (bit & 31)) & 1u : 0u;
              unsigned long long c = g;
              for(int t = 0; t <= L && c; ++t)
                {
                  unsigned long long v = (unsigned long long)r[t] + c;
                  r[t] = uint32_t(v);
                  c = v >> 32;
                }
            }
          if((r[L - 1] >> (NB - 32 * (L - 1))) == 0 && r[L] == 0)
            break;
          ++sh;
        }
      if(neg)
        {
          unsigned long long c = 1;
          for(int t = 0; t < L; ++t)
            {
              unsigned long long v = (unsigned long long)(~r[t]) + c;
              r[t] = uint32_t(v);
              c = v >> 32;
            }
        }
      for(int t = 0; t < L; ++t)
        om[(int64_t)t * total + eo] = int32_t(r[t]);
      oe[eo] = ebase + sh;
    }

    __global__ void normalize_blocks(const int32_t *__restrict__ Ci, int W,
                                     int64_t in_plane, int64_t mp, int64_t m,
                                     int64_t n, int64_t total,
                                     const int64_t *__restrict__ EA,
                                     const int64_t *__restrict__ EB, int L,
                                     int32_t *__restrict__ om,
                                     int64_t *__restrict__ oe, int64_t j0 = 0,
                                     int64_t pc = -1)
    {
      if(pc < 0)
        pc = n;
      const int64_t q = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      const int64_t batch = total / (m * n);
      if(q >= batch * m * pc)
        return;
      const int64_t b = q / (m * pc), rem = q % (m * pc);
      const int64_t jj = rem / m, i = rem % m, j = j0 + jj;
      round_entry(Ci, W, in_plane, b * mp * pc + jj * mp + i, L, total,
                  b * m * n + j * m + i, EA[b * m + i] + EB[b * n + j], om, oe);
    }

    // the same rounding with the limbs in registers: static indices, a
    // barrel shifter over limbs and funnel shifts instead of local arrays
    // with data-dependent indices; the rare cases (fewer bits than a float
    // holds, or a rounding carry into the top bit) take round_entry
    template <int WM, int LT>
    __global__ void normalize_blocks_reg(const int32_t *__restrict__ Ci, int W,
                                         int64_t in_plane, int64_t mp, int64_t m,
                                         int64_t n, int64_t total,
                                         const int64_t *__restrict__ EA,
                                         const int64_t *__restrict__ EB,
                                         int32_t *__restrict__ om,
                                         int64_t *__restrict__ oe, int64_t j0,
                                         int64_t pc)
    {
      const int64_t q = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      const int64_t batch = total / (m * n);
      if(q >= batch * m * pc)
        return;
      const int64_t b = q / (m * pc), rem = q % (m * pc);
      const int64_t jj = rem / m, i = rem % m, j = j0 + jj;
      const int64_t ei = b * mp * pc + jj * mp + i;
      const int64_t eo = b * m * n + j * m + i;
      const int64_t ebase = EA[b * m + i] + EB[b * n + j];
      uint32_t a[WM];
#pragma unroll
      for(int t = 0; t < WM; ++t)
        a[t] = t < W ? uint32_t(Ci[(int64_t)t * in_plane + ei]) : 0u;
      const bool neg = int32_t(Ci[(int64_t)(W - 1) * in_plane + ei]) < 0;
      if(neg)
        {
          unsigned long long c = 1;
#pragma unroll
          for(int t = 0; t < WM; ++t)
            if(t < W)
              {
                const unsigned long long v = (unsigned long long)(~a[t]) + c;
                a[t] = uint32_t(v);
                c = v >> 32;
              }
        }
      int bl = 0;
#pragma unroll
      for(int t = 0; t < WM; ++t)
        if(a[t])
          bl = 32 * t + 32 - __clz(a[t]);
      if(bl == 0)
        {
#pragma unroll
          for(int t = 0; t < LT; ++t)
            om[(int64_t)t * total + eo] = 0;
          oe[eo] = MPF_EZ;
          return;
        }
      constexpr int NB = 32 * LT - 2;
      const int sh = bl - NB;
      if(sh < 1)
        {
          round_entry(Ci, W, in_plane, ei, LT, total, eo, ebase, om, oe);
          return;
        }
      // the rounding bit, bit sh - 1
      const int gq = (sh - 1) >> 5, gr = (sh - 1) & 31;
      uint32_t g = 0;
#pragma unroll
      for(int t = 0; t < WM; ++t)
        if(t == gq)
          g = (a[t] >> gr) & 1u;
      // a >>= sh: whole limbs by a barrel shifter, then the bits
      const int qs = sh >> 5, rr = sh & 31;
#pragma unroll
      for(int bb = 0; (1 << bb) < WM; ++bb)
        if(qs & (1 << bb))
          {
#pragma unroll
            for(int t = 0; t < WM; ++t)
              a[t] = (t + (1 << bb) < WM) ? a[t + (1 << bb)] : 0u;
          }
      uint32_t r[LT + 1];
#pragma unroll
      for(int t = 0; t <= LT; ++t)
        r[t] = __funnelshift_r(a[t], (t + 1 < WM) ? a[t + 1] : 0u, rr);
      unsigned long long c = g;
#pragma unroll
      for(int t = 0; t <= LT; ++t)
        {
          const unsigned long long v = (unsigned long long)r[t] + c;
          r[t] = uint32_t(v);
          c = v >> 32;
        }
      if((r[LT - 1] >> (NB - 32 * (LT - 1))) != 0 || r[LT] != 0)
        {
          round_entry(Ci, W, in_plane, ei, LT, total, eo, ebase, om, oe);
          return;
        }
      if(neg)
        {
          unsigned long long cc = 1;
#pragma unroll
          for(int t = 0; t < LT; ++t)
            {
              const unsigned long long v = (unsigned long long)(~r[t]) + cc;
              r[t] = uint32_t(v);
              cc = v >> 32;
            }
        }
#pragma unroll
      for(int t = 0; t < LT; ++t)
        om[(int64_t)t * total + eo] = int32_t(r[t]);
      oe[eo] = ebase + sh;
    }
  }

  Mpf_Gram::Mpf_Gram(int N, int L, int digits, int64_t rows_total,
                     Gram_Options options, cudaStream_t stream)
      : N_(N), L_(L), D_(digits),
        G_(N, digits, rows_total + 16 * 4096, options, stream),
        stream_(stream), timing_(options.timing)
  {
    // (the constructor of G_ received rows_total + 16 * 4096: each chunk's
    // rows are padded to a multiple of 16)
    if(L > LMAX)
      throw std::invalid_argument("Mpf_Gram: too many limbs");
    // the column maximum is shifted by at most 32 bits to the left
    // (s = e_ref - e >= 32L + 2 - 8 digits >= -32)
    if(8 * digits > 32 * L + 34)
      throw std::invalid_argument("Mpf_Gram: digits > (32 L + 34) / 8");
    if(G_.limbs_out() > 96)
      throw std::invalid_argument("Mpf_Gram: too many output limbs");
    // residue tables: c[k] = 2^(32k) mod m_i (k = 0..L+2), inv2[r] = 2^-r mod m_i
    W_ = L + 3;
    const auto &pr = G_.basis().primes;
    std::vector<uint32_t> ct(pr.size() * W_), cp(pr.size() * 32);
    for(size_t q = 0; q < pr.size(); ++q)
      {
        unsigned long long v = 1, m = pr[q], h = (m + 1) / 2, x = 1;
        for(int w = 0; w < W_; ++w)
          {
            ct[q * W_ + w] = uint32_t(v);
            v = (v << 32) % m;
          }
        for(int r = 0; r < 32; ++r)
          {
            cp[q * 32 + r] = uint32_t(x);
            x = x * h % m; // h = 2^-1 mod m (m odd)
          }
      }
    MG_CHECK(cudaMalloc(&d_ctab_, ct.size() * 4));
    MG_CHECK(cudaMalloc(&d_ctop_, cp.size() * 4));
    MG_CHECK(cudaMemcpy(d_ctab_, ct.data(), ct.size() * 4, cudaMemcpyHostToDevice));
    MG_CHECK(cudaMemcpy(d_ctop_, cp.data(), cp.size() * 4, cudaMemcpyHostToDevice));
    MG_CHECK(cudaMalloc(&d_top_, N * sizeof(int64_t)));
    MG_CHECK(cudaMalloc(&d_Eref_, N * sizeof(int64_t)));
    MG_CHECK(cudaMemsetAsync(d_top_, 0, N * sizeof(int64_t), stream_));
    MG_CHECK(cudaMalloc(&d_Qint_, size_t(G_.limbs_out()) * G_.packed_size() * 4));
    MG_CHECK(cudaEventCreate(&ev0_));
    MG_CHECK(cudaEventCreate(&ev1_));
  }

  Mpf_Gram::~Mpf_Gram()
  {
    cudaFree(d_top_);
    cudaFree(d_Eref_);
    cudaFree(d_ctab_);
    cudaFree(d_ctop_);
    cudaFree(d_Qint_);
    cudaEventDestroy(ev0_);
    cudaEventDestroy(ev1_);
  }

  void Mpf_Gram::tic()
  {
    if(timing_)
      MG_CHECK(cudaEventRecord(ev0_, stream_));
  }

  void Mpf_Gram::toc(const std::string &phase)
  {
    if(!timing_)
      return;
    MG_CHECK(cudaEventRecord(ev1_, stream_));
    MG_CHECK(cudaEventSynchronize(ev1_));
    float ms = 0;
    MG_CHECK(cudaEventElapsedTime(&ms, ev0_, ev1_));
    ms_[phase] += ms;
  }

  void Mpf_Gram::scan(const Mpf_Matrix &P)
  {
    if(P.cols != N_ || P.L != L_)
      throw std::invalid_argument("Mpf_Gram::scan: shape");
    tic();
    col_tops<<<dim3(unsigned(P.cols), unsigned((P.rows + 4095) / 4096)), 256,
               0, stream_>>>(
      P.m, P.e, L_, P.rows, P.cols, P.rows * P.cols,
      reinterpret_cast<unsigned long long *>(d_top_));
    MG_CHECK(cudaGetLastError());
    refs_ready_ = false;
    toc("column maxima");
  }

  void Mpf_Gram::add(const Mpf_Matrix &P)
  {
    if(P.cols != N_ || P.L != L_)
      throw std::invalid_argument("Mpf_Gram::add: shape");
    if(!refs_ready_)
      {
        make_refs<<<nblocks(N_), threads, 0, stream_>>>(
          reinterpret_cast<unsigned long long *>(d_top_), d_Eref_, N_, D_);
        MG_CHECK(cudaGetLastError());
        refs_ready_ = true;
      }
    const int64_t ld = (P.rows + 15) / 16 * 16;
    G_.add_residues(ld, [&](int i0, int gc, int8_t *R, int64_t plane) {
      const size_t sm = size_t(gc) * W_ * 4;
#define MPF_RES(LT)                                                           \
  mpf_residues<LT><<<nblocks(ld * N_), threads, sm, stream_>>>(              \
    P.m, P.e, L_, P.rows, P.cols, ld, P.rows * P.cols, 0, d_Eref_,            \
    d_ctab_ + int64_t(i0) * W_,                                               \
    d_ctop_ + int64_t(i0) * 32, G_.device_primes() + i0,                      \
    G_.device_minv() + i0, gc, R, plane)
      if(L_ == 28)
        MPF_RES(28);
      else if(L_ == 30)
        MPF_RES(30);
      else
        MPF_RES(0);
#undef MPF_RES
      MG_CHECK(cudaGetLastError());
    });
  }

  void Mpf_Gram::finish(int32_t *d_m, int64_t *d_e)
  {
    G_.finish(d_Qint_);
    tic();
    normalize_packed<<<nblocks(G_.packed_size()), threads, 0, stream_>>>(
      d_Qint_, G_.limbs_out(), G_.packed_size(), N_, d_Eref_, L_, d_m, d_e);
    MG_CHECK(cudaGetLastError());
    MG_CHECK(cudaMemsetAsync(d_top_, 0, N_ * sizeof(int64_t), stream_));
    refs_ready_ = false;
    toc("round Q to floats");
  }

  Mpf_Gemm::Mpf_Gemm(int m, int n, int k, int batch, int L, int digits,
                     Gram_Options options, int panel, cudaStream_t stream)
      : m_(m), n_(n), k_(k), kp_((k + 15) / 16 * 16), mp_((m + 3) / 4 * 4),
        batch_(batch), L_(L), D_(digits), fast_rounding_(options.fast_rounding),
        G_(mp_, n, kp_, batch, digits, digits, options, panel, stream),
        stream_(stream), timing_(options.timing)
  {
    if(L > LMAX || 8 * digits > 32 * L + 34 || G_.limbs_out() > 96)
      throw std::invalid_argument("Mpf_Gemm: sizes");
    const auto &pr = G_.basis().primes;
    const int CK = L + 3;
    std::vector<uint32_t> ct(pr.size() * CK);
    for(size_t q = 0; q < pr.size(); ++q)
      {
        unsigned long long v = 1, mm = pr[q];
        for(int w = 0; w < CK; ++w)
          {
            ct[q * CK + w] = uint32_t(v);
            v = (v << 32) % mm;
          }
      }
    MG_CHECK(cudaMalloc(&d_ctab_, ct.size() * 4));
    MG_CHECK(cudaMemcpy(d_ctab_, ct.data(), ct.size() * 4, cudaMemcpyHostToDevice));
    MG_CHECK(cudaMalloc(&d_topA_, size_t(batch) * m * 8));
    MG_CHECK(cudaMalloc(&d_topB_, size_t(batch) * n * 8));
    MG_CHECK(cudaMalloc(&d_EA_, size_t(batch) * m * 8));
    MG_CHECK(cudaMalloc(&d_EB_, size_t(batch) * n * 8));
    MG_CHECK(cudaEventCreate(&ev0_));
    MG_CHECK(cudaEventCreate(&ev1_));
  }

  Mpf_Gemm::~Mpf_Gemm()
  {
    for(void *p : {(void *)d_ctab_, (void *)d_topA_, (void *)d_topB_,
                   (void *)d_EA_, (void *)d_EB_, (void *)d_exact_})
      cudaFree(p);
    cudaEventDestroy(ev0_);
    cudaEventDestroy(ev1_);
  }

  void Mpf_Gemm::tic()
  {
    if(timing_)
      MG_CHECK(cudaEventRecord(ev0_, stream_));
  }

  void Mpf_Gemm::toc(const std::string &phase)
  {
    if(!timing_)
      return;
    MG_CHECK(cudaEventRecord(ev1_, stream_));
    MG_CHECK(cudaEventSynchronize(ev1_));
    float ms = 0;
    MG_CHECK(cudaEventElapsedTime(&ms, ev0_, ev1_));
    ms_[phase] += ms;
  }

  void Mpf_Gemm::multiply(const Mpf_Matrix &A, const Mpf_Matrix &B,
                          int32_t *d_cm, int64_t *d_ce, int triangles)
  {
    if(A.rows != k_ || A.cols != m_ || B.rows != k_ || B.cols != n_
       || A.L != L_ || B.L != L_)
      throw std::invalid_argument("Mpf_Gemm::multiply: shapes");
    tic();
    MG_CHECK(cudaMemsetAsync(d_topA_, 0, size_t(batch_) * m_ * 8, stream_));
    MG_CHECK(cudaMemsetAsync(d_topB_, 0, size_t(batch_) * n_ * 8, stream_));
    const int64_t plA = int64_t(batch_) * k_ * m_, plB = int64_t(batch_) * k_ * n_;
    col_tops<<<dim3(unsigned(m_), unsigned((k_ + 4095) / 4096), unsigned(batch_)),
               256, 0, stream_>>>(A.m, A.e, L_, k_, m_, plA,
                                  reinterpret_cast<unsigned long long *>(d_topA_),
                                  int64_t(k_) * m_, int64_t(m_));
    col_tops<<<dim3(unsigned(n_), unsigned((k_ + 4095) / 4096), unsigned(batch_)),
               256, 0, stream_>>>(B.m, B.e, L_, k_, n_, plB,
                                  reinterpret_cast<unsigned long long *>(d_topB_),
                                  int64_t(k_) * n_, int64_t(n_));
    make_refs<<<nblocks(int64_t(batch_) * m_), threads, 0, stream_>>>(
      reinterpret_cast<unsigned long long *>(d_topA_), d_EA_, int64_t(batch_) * m_, D_);
    make_refs<<<nblocks(int64_t(batch_) * n_), threads, 0, stream_>>>(
      reinterpret_cast<unsigned long long *>(d_topB_), d_EB_, int64_t(batch_) * n_, D_);
    MG_CHECK(cudaGetLastError());
    toc("column maxima");
    const int CK = L_ + 3;
    auto launch = [&](const int32_t *mm, const int64_t *ee, int64_t cols,
                      int64_t pl, int64_t in_bs, const int64_t *Eref,
                      int64_t eref_bs, int64_t r_bs, int i0, int gc, int8_t *R,
                      int64_t plane) {
      const size_t sm = size_t(gc) * CK * 4;
      const int64_t nt = int64_t(kp_) * cols;
      const dim3 grid(unsigned(nblocks(nt)), unsigned(batch_));
#define MPF_RES2(LT)                                                          \
  mpf_residues<LT><<<grid, threads, sm, stream_>>>(                          \
    mm, ee, L_, k_, cols, kp_, pl, 0, Eref, d_ctab_ + int64_t(i0) * CK,       \
    nullptr, G_.device_primes() + i0, G_.device_minv() + i0, gc, R, plane,    \
    in_bs, eref_bs, r_bs)
      if(L_ == 28)
        MPF_RES2(28);
      else if(L_ == 30)
        MPF_RES2(30);
      else
        MPF_RES2(0);
#undef MPF_RES2
      MG_CHECK(cudaGetLastError());
    };
    G_.multiply_with(
      [&](int i0, int gc, int8_t *RA, int64_t planeA) {
        if(mp_ != m_)
          MG_CHECK(cudaMemsetAsync(RA, 0, size_t(3 * gc) * planeA, stream_));
        launch(A.m, A.e, m_, plA, int64_t(k_) * m_, d_EA_, int64_t(m_),
               int64_t(kp_) * mp_, i0, gc, RA, planeA);
      },
      [&](int i0, int gc, int j0, int pc, int8_t *RB, int64_t planeB) {
        launch(B.m + int64_t(j0) * k_, B.e + int64_t(j0) * k_, pc, plB,
               int64_t(k_) * n_, d_EB_ + j0, int64_t(n_), int64_t(kp_) * pc, i0,
               gc, RB, planeB);
      },
      nullptr,
      [&](const int32_t *panel, int64_t cplane, int j0, int pc) {
        const int64_t total = int64_t(batch_) * m_ * n_;
        if(fast_rounding_ && L_ == 28 && G_.limbs_out() <= 64)
          normalize_blocks_reg<64, 28><<<nblocks(int64_t(batch_) * m_ * pc),
                                         threads, 0, stream_>>>(
            panel, G_.limbs_out(), cplane, mp_, m_, n_, total, d_EA_, d_EB_,
            d_cm, d_ce, j0, pc);
        else
          normalize_blocks<<<nblocks(int64_t(batch_) * m_ * pc), threads, 0,
                             stream_>>>(panel, G_.limbs_out(), cplane, mp_, m_,
                                        n_, total, d_EA_, d_EB_, L_, d_cm, d_ce,
                                        j0, pc);
        MG_CHECK(cudaGetLastError());
      }, triangles);
  }
}
