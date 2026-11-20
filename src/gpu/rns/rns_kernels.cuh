// GPU backend for SDPB: device kernels and the cuBLASLt INT8 product shared
// by the residue-arithmetic classes (Gram_Accumulator, Rns_Gemm).
//
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: sdpvibe/rns.py (residues as an INT8 product
// with a table of 256^j mod m_i; CRT as an INT8 product with the digits of
// M/m_i, then carries).  Portions adapted from SDPviBe, Copyright (c) 2026
// Antoine Tilloy, MIT License (docs/GPU_backend_credits.md).
#pragma once

#include <cublasLt.h>
#include <cuda_runtime.h>

#include <array>
#include <cstdint>
#include <map>
#include <sstream>
#include <stdexcept>

#define RNS_CUDA_CHECK(x)                                                     \
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
  // cudaMallocAsync that, when the device is out of memory, waits for the
  // stream, returns the stream-ordered pool's unused memory to the device and
  // tries once more (memory freed with cudaFreeAsync stays in the pool until
  // a synchronization)
  inline cudaError_t malloc_async_retry(void **p, size_t bytes, cudaStream_t s)
  {
    cudaError_t e = cudaMallocAsync(p, bytes, s);
    if(e != cudaErrorMemoryAllocation)
      return e;
    (void)cudaGetLastError();
    cudaStreamSynchronize(s);
    int dev = 0;
    cudaGetDevice(&dev);
    cudaMemPool_t pool;
    if(cudaDeviceGetDefaultMemPool(&pool, dev) == cudaSuccess)
      cudaMemPoolTrimTo(pool, 0);
    return cudaMallocAsync(p, bytes, s);
  }
}

#define RNS_LT_CHECK(x)                                                       \
  do                                                                          \
    {                                                                         \
      cublasStatus_t s_ = (x);                                                \
      if(s_ != CUBLAS_STATUS_SUCCESS)                                         \
        {                                                                     \
          std::ostringstream ss_;                                             \
          ss_ << __FILE__ << ":" << __LINE__ << ": cuBLASLt status "         \
              << int(s_);                                                     \
          throw std::runtime_error(ss_.str());                                \
        }                                                                     \
  } while(0)

namespace sdpb_gpu
{
  namespace rns_detail
  {
  namespace // internal linkage: each translation unit gets its own kernels
  {
    constexpr int threads = 256;
    inline int64_t blocks_for(int64_t n) { return (n + threads - 1) / threads; }

    // v mod m in [0, m) for |v| < 2^50: float64 quotient, one correction
    __device__ __forceinline__ int32_t mod_d(long long v, int32_t m,
                                             double minv)
    {
      long long q = (long long)floor((double)v * minv);
      long long r = v - q * m;
      if(r < 0)
        r += m;
      else if(r >= m)
        r -= m;
      return (int32_t)r;
    }

    // C (m_ent x 2g, column-major): column 2p = sum_t d_t T_hi(p, t), 2p+1 the
    // lo part.  r = 128 z_hi + z_lo mod m_p, centred, r = 128 h + l with
    // l in [-64, 63], h in [-64, 64], s = h + l in [-127, 127]; planes
    // R[(X*g + p)*plane + e_off + e] for X = h, l, s.
    __global__ void residue_finish(const int32_t *__restrict__ C,
                                          int64_t m_ent, int g,
                                          const int32_t *__restrict__ primes,
                                          const double *__restrict__ minv,
                                          int8_t *__restrict__ R,
                                          int64_t plane, int64_t e_off)
    {
      int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(e >= m_ent)
        return;
      for(int p = 0; p < g; ++p)
        {
          long long z
            = 128LL * C[(2 * p) * m_ent + e] + C[(2 * p + 1) * m_ent + e];
          int32_t m = primes[p];
          int32_t r = mod_d(z, m, minv[p]);
          if(r > m / 2)
            r -= m;
          int32_t lo = ((r + 64) & 127) - 64;
          int32_t hi = (r - lo) >> 7;
          R[(int64_t)(0 * g + p) * plane + e_off + e] = (int8_t)hi;
          R[(int64_t)(1 * g + p) * plane + e_off + e] = (int8_t)lo;
          R[(int64_t)(2 * g + p) * plane + e_off + e] = (int8_t)(hi + lo);
        }
    }

    // CRT inputs: y_i = q_i (M/m_i)^-1 mod m_i, centred, y = 128 yh + yl;
    // k = round(sum y_i / m_i).  Qres[i*qplane + e0 + e]; Y rows of rpad bytes.
    __global__ void crt_pre(const uint16_t *__restrict__ Qres,
                                   int64_t qplane, int64_t e0, int64_t ne,
                                   int64_t ne_pad, int r, int rpad,
                                   const int32_t *__restrict__ primes,
                                   const double *__restrict__ minv,
                                   const int32_t *__restrict__ inv,
                                   int8_t *__restrict__ Yh,
                                   int8_t *__restrict__ Yl,
                                   int32_t *__restrict__ K, int ystride = 0)
    {
      // ystride: bytes between entries' rows of Yh and Yl (default rpad)
      if(ystride == 0)
        ystride = rpad;
      int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(e >= ne_pad)
        return;
      double ks = 0;
      for(int i0 = 0; i0 < rpad; i0 += 16)
        {
          alignas(16) int8_t bh[16];
          alignas(16) int8_t bl[16];
          for(int u = 0; u < 16; ++u)
            {
              int i = i0 + u;
              int32_t yh = 0, yl = 0;
              if(i < r && e < ne)
                {
                  int32_t m = primes[i];
                  long long q = Qres[(int64_t)i * qplane + e0 + e];
                  int32_t y = mod_d(q * inv[i], m, minv[i]);
                  if(y > m / 2)
                    y -= m;
                  ks += (double)y * minv[i];
                  yl = ((y + 64) & 127) - 64;
                  yh = (y - yl) >> 7;
                }
              bh[u] = (int8_t)yh;
              bl[u] = (int8_t)yl;
            }
          *reinterpret_cast<int4 *>(Yh + e * ystride + i0)
            = *reinterpret_cast<const int4 *>(bh);
          *reinterpret_cast<int4 *>(Yl + e * ystride + i0)
            = *reinterpret_cast<const int4 *>(bl);
        }
      K[e] = (int32_t)llrint(ks);
    }

    // x = sum_j (Z_j - k Mdig_j) 128^j (base-128 single-product CRT,
    // Rns_Basis::crt_tables_128), carried into 7-bit digits and packed into W
    // two's complement 32-bit limbs out[w*out_plane + e0 + e]; the digits past
    // E2 are the sign
    __global__ void crt_carry_128(const int32_t *__restrict__ Z,
                                  const int32_t *__restrict__ K,
                                  const int32_t *__restrict__ Mdig, int64_t ne,
                                  int64_t ldz, int E2, int W,
                                  int32_t *__restrict__ out, int64_t out_plane,
                                  int64_t e0, int *err,
                                  const float *__restrict__ Kf = nullptr)
    {
      int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(e >= ne)
        return;
      long long v = 0, k = Kf ? (long long)__float2int_rn(Kf[e]) : K[e];
      // 32 digits of 7 bits make exactly 7 words: digits in groups of 32,
      // each at a static bit position
      const int nd = max(E2, (32 * W + 6) / 7), ng = (nd + 31) / 32;
      bool bad = false;
      uint32_t sign = 0;
      for(int gq = 0; gq < ng; ++gq)
        {
          uint32_t wd[7] = {0, 0, 0, 0, 0, 0, 0};
#pragma unroll
          for(int u = 0; u < 32; ++u)
            {
              const int j = gq * 32 + u;
              if(j < E2)
                v += (long long)Z[(int64_t)j * ldz + e] - k * Mdig[j];
              const uint32_t dig = uint32_t(v & 127);
              v = (v - (long long)dig) >> 7;
              const int bit = 7 * u, wq = bit >> 5, off = bit & 31;
              wd[wq] |= dig << off;
              if(off > 25)
                wd[wq + 1] |= dig >> (32 - off);
            }
#pragma unroll
          for(int q = 0; q < 7; ++q)
            {
              const int wi = gq * 7 + q;
              if(wi < W)
                {
                  out[(int64_t)wi * out_plane + e0 + e] = (int32_t)wd[q];
                  if(wi == W - 1)
                    sign = (wd[q] >> 31) ? 0xffffffffu : 0u;
                }
              else
                bad = bad || wd[q] != sign;
            }
        }
      if(bad || !((v == 0 && sign == 0) || (v == -1 && sign != 0)))
        atomicExch(err, 1);
    }

    // x = sum_j (128 Zh_j + Zl_j - k Mdig_j) 256^j, carried into bytes,
    // written as two's complement 32-bit limbs out[w*out_plane + e0 + e]
    __global__ void crt_carry(const int32_t *__restrict__ Zh,
                                     const int32_t *__restrict__ Zl,
                                     const int32_t *__restrict__ K,
                                     const int32_t *__restrict__ Mdig,
                                     int64_t ne, int64_t ldz, int E,
                                     int32_t *__restrict__ out,
                                     int64_t out_plane, int64_t e0, int *err,
                                     const float *__restrict__ Kf = nullptr)
    {
      int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(e >= ne)
        return;
      long long v = 0, k = Kf ? (long long)__float2int_rn(Kf[e]) : K[e];
      uint32_t word = 0;
      for(int j = 0; j < E; ++j)
        {
          v += (Zl ? 128LL * Zh[(int64_t)j * ldz + e] + Zl[(int64_t)j * ldz + e]
                   : (long long)Zh[(int64_t)j * ldz + e])
               - k * Mdig[j];
          uint32_t byte = (uint32_t)(v & 255);
          v = (v - (long long)byte) >> 8;
          word |= byte << (8 * (j & 3));
          if((j & 3) == 3)
            {
              out[(int64_t)(j >> 2) * out_plane + e0 + e] = (int32_t)word;
              word = 0;
            }
        }
      bool neg
        = (int32_t)(out[(int64_t)(E / 4 - 1) * out_plane + e0 + e]) < 0;
      if(!((v == 0 && !neg) || (v == -1 && neg)))
        atomicExch(err, 1);
    }
  }
  }

  // C (m x n, int32) = A^T B with A (k x m), B (k x n) int8, column-major
  // ("TN", the layout of cuBLASLt's INT8 tensor-core kernels), optionally
  // strided-batched.  Returns the INT8 operations issued (2 per MAC).
  inline double lt_gemm_tn(cublasLtHandle_t lt, void *workspace,
                           size_t workspace_bytes, cudaStream_t stream, int m,
                           int n, int k, const int8_t *A, int lda,
                           int64_t strideA, const int8_t *B, int ldb,
                           int64_t strideB, int32_t *C, int ldc,
                           int64_t strideC, int batch, int32_t alpha = 1,
                           int32_t beta = 0)
  {
    cublasLtMatmulDesc_t op;
    RNS_LT_CHECK(cublasLtMatmulDescCreate(&op, CUBLAS_COMPUTE_32I, CUDA_R_32I));
    cublasOperation_t tA = CUBLAS_OP_T, tB = CUBLAS_OP_N;
    RNS_LT_CHECK(cublasLtMatmulDescSetAttribute(
      op, CUBLASLT_MATMUL_DESC_TRANSA, &tA, sizeof(tA)));
    RNS_LT_CHECK(cublasLtMatmulDescSetAttribute(
      op, CUBLASLT_MATMUL_DESC_TRANSB, &tB, sizeof(tB)));
    cublasLtMatrixLayout_t La, Lb, Lc;
    RNS_LT_CHECK(cublasLtMatrixLayoutCreate(&La, CUDA_R_8I, k, m, lda));
    RNS_LT_CHECK(cublasLtMatrixLayoutCreate(&Lb, CUDA_R_8I, k, n, ldb));
    RNS_LT_CHECK(cublasLtMatrixLayoutCreate(&Lc, CUDA_R_32I, m, n, ldc));
    if(batch > 1)
      {
        cublasLtMatrixLayout_t Ls[3] = {La, Lb, Lc};
        int64_t ss[3] = {strideA, strideB, strideC};
        for(int q = 0; q < 3; ++q)
          {
            RNS_LT_CHECK(cublasLtMatrixLayoutSetAttribute(
              Ls[q], CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batch,
              sizeof(batch)));
            RNS_LT_CHECK(cublasLtMatrixLayoutSetAttribute(
              Ls[q], CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET, &ss[q],
              sizeof(ss[q])));
          }
      }
    // heuristic results cached by shape (the query costs ~0.1 ms per call)
    static std::map<std::array<int64_t, 11>, cublasLtMatmulHeuristicResult_t> cache;
    const std::array<int64_t, 11> key{m, n, k, lda, ldb, ldc, strideA,
                                      strideB, strideC, batch,
                                      int64_t(workspace_bytes)};
    cublasLtMatmulHeuristicResult_t res{};
    int nres = 0;
    cublasLtMatmulPreference_t pref = nullptr;
    auto hit = cache.find(key);
    if(hit != cache.end())
      {
        res = hit->second;
        nres = 1;
      }
    else
      {
        RNS_LT_CHECK(cublasLtMatmulPreferenceCreate(&pref));
        RNS_LT_CHECK(cublasLtMatmulPreferenceSetAttribute(
          pref, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &workspace_bytes,
          sizeof(workspace_bytes)));
        RNS_LT_CHECK(cublasLtMatmulAlgoGetHeuristic(lt, op, La, Lb, Lc, Lc,
                                                    pref, 1, &res, &nres));
        if(nres > 0)
          cache[key] = res;
      }
    if(nres == 0)
      {
        std::ostringstream ss;
        ss << "cuBLASLt: no INT8 algorithm for m=" << m << " n=" << n
           << " k=" << k << " lda=" << lda << " ldb=" << ldb
           << " ldc=" << ldc << " batch=" << batch;
        throw std::runtime_error(ss.str());
      }
    RNS_LT_CHECK(cublasLtMatmul(lt, op, &alpha, A, La, B, Lb, &beta, C, Lc,
                                C, Lc, &res.algo, workspace, workspace_bytes,
                                stream));
    if(pref)
      cublasLtMatmulPreferenceDestroy(pref);
    cublasLtMatrixLayoutDestroy(La);
    cublasLtMatrixLayoutDestroy(Lb);
    cublasLtMatrixLayoutDestroy(Lc);
    cublasLtMatmulDescDestroy(op);
    return 2.0 * m * n * double(k) * batch;
  }
}

#include <mma.h>

#include <algorithm>

namespace sdpb_gpu
{
  namespace rns_detail
  {
  namespace
  {
    // Fused residues (no int32 intermediate): for 16-entry tiles, the INT8
    // tensor-core product digits (16 x NB) @ T (NB x 2G) runs in WMMA
    // fragments, and the modular finish happens in shared memory; h, l, s are
    // written as int8 planes R[(X*G + p)*plane + e_off + e].  G <= 32.
    // Digits: entry e at d + e*NB (NB % 4 == 0); T: [G][2][NB] (hi, lo rows).
    template <int G>
    __global__ void residues_fused(const int8_t *__restrict__ d,
                                   int64_t entries, int NB,
                                   const int8_t *__restrict__ T,
                                   const int32_t *__restrict__ primes,
                                   const double *__restrict__ minv,
                                   int8_t *__restrict__ R, int64_t plane,
                                   int64_t e_off, int gc)
    {
      using namespace nvcuda;
      constexpr int WARPS = 8, TILE = 16, NC = 2 * G;
      const int KP = (NB + 15) / 16 * 16;
      extern __shared__ __align__(32) unsigned char smem[];
      // k split in 16-byte chunks so that every fragment starts 32-byte aligned
      int8_t *sT = reinterpret_cast<int8_t *>(smem);                   // [KP/16][NC][16]
      int8_t *sD = sT + NC * KP;                                        // [WARPS][KP/16][16][16]
      int32_t *sZ = reinterpret_cast<int32_t *>(sD + WARPS * TILE * KP); // [WARPS][16][NC]
      const int warp = threadIdx.x / 32, lane = threadIdx.x % 32;
      // table, zero padded in k and in columns >= 2 gc
      for(int q = threadIdx.x; q < NC * KP; q += blockDim.x)
        {
          int c = q / KP, t = q % KP;
          sT[((t / 16) * NC + c) * 16 + t % 16]
            = (c < 2 * gc && t < NB) ? T[c * NB + t] : 0;
        }
      __syncthreads();
      const int64_t tiles = (entries + TILE - 1) / TILE;
      for(int64_t tile = int64_t(blockIdx.x) * WARPS + warp; tile < tiles;
          tile += int64_t(gridDim.x) * WARPS)
        {
          const int64_t e0 = tile * TILE;
          int8_t *myD = sD + warp * TILE * KP;
          // stage 16 entries' digits (4-byte words), zero padded
          for(int q = lane; q < TILE * (KP / 4); q += 32)
            {
              int row = q / (KP / 4), w = q % (KP / 4);
              int32_t v = 0;
              if(e0 + row < entries && 4 * w < NB)
                v = *reinterpret_cast<const int32_t *>(d + (e0 + row) * NB + 4 * w);
              int t = 4 * w;
              *reinterpret_cast<int32_t *>(myD + ((t / 16) * TILE + row) * 16 + t % 16) = v;
            }
          __syncwarp();
          wmma::fragment<wmma::accumulator, 16, 16, 16, int> acc[NC / 16];
          for(int c = 0; c < NC / 16; ++c)
            wmma::fill_fragment(acc[c], 0);
          for(int k0 = 0; k0 < KP; k0 += 16)
            {
              wmma::fragment<wmma::matrix_a, 16, 16, 16, signed char,
                             wmma::row_major>
                a;
              wmma::load_matrix_sync(a, myD + (k0 / 16) * 256, 16);
              for(int c = 0; c < NC / 16; ++c)
                {
                  wmma::fragment<wmma::matrix_b, 16, 16, 16, signed char,
                                 wmma::col_major>
                    b;
                  wmma::load_matrix_sync(b, sT + ((k0 / 16) * NC + c * 16) * 16, 16);
                  wmma::mma_sync(acc[c], a, b, acc[c]);
                }
            }
          int32_t *myZ = sZ + warp * TILE * NC;
          for(int c = 0; c < NC / 16; ++c)
            wmma::store_matrix_sync(myZ + c * 16, acc[c], NC,
                                    wmma::mem_row_major);
          __syncwarp();
          // finish: lane -> (prime p, entries), consecutive lanes on entries
          for(int q = lane; q < TILE * gc; q += 32)
            {
              int row = q % TILE, p = q / TILE;
              int64_t e = e0 + row;
              if(e >= entries)
                continue;
              long long z = 128LL * myZ[row * NC + 2 * p] + myZ[row * NC + 2 * p + 1];
              int32_t m = primes[p];
              int32_t r = mod_d(z, m, minv[p]);
              if(r > m / 2)
                r -= m;
              int32_t lo = ((r + 64) & 127) - 64;
              int32_t hi = (r - lo) >> 7;
              R[(int64_t)(0 * gc + p) * plane + e_off + e] = (int8_t)hi;
              R[(int64_t)(1 * gc + p) * plane + e_off + e] = (int8_t)lo;
              R[(int64_t)(2 * gc + p) * plane + e_off + e] = (int8_t)(hi + lo);
            }
          __syncwarp();
        }
    }

    // launch helper: smem = NC*KP + 8*16*KP + 8*16*NC*4 bytes
    inline void launch_residues_fused(const int8_t *d, int64_t entries, int NB,
                                      const int8_t *T, const int32_t *primes,
                                      const double *minv, int8_t *R,
                                      int64_t plane, int64_t e_off, int gc,
                                      cudaStream_t stream)
    {
      const int KP = (NB + 15) / 16 * 16;
      const int64_t tiles = (entries + 15) / 16;
      int grid = int(std::min<int64_t>((tiles + 7) / 8, 132 * 8));
      if(gc <= 8)
        {
          size_t sm = 16 * KP + 8 * 16 * KP + 8 * 16 * 16 * 4;
          residues_fused<8><<<grid, 256, sm, stream>>>(d, entries, NB, T, primes, minv, R, plane, e_off, gc);
        }
      else if(gc <= 16)
        {
          size_t sm = 32 * KP + 8 * 16 * KP + 8 * 16 * 32 * 4;
          residues_fused<16><<<grid, 256, sm, stream>>>(d, entries, NB, T, primes, minv, R, plane, e_off, gc);
        }
      else
        {
          size_t sm = 64 * KP + 8 * 16 * KP + 8 * 16 * 64 * 4;
          if(sm > 48 * 1024)
            RNS_CUDA_CHECK(cudaFuncSetAttribute(residues_fused<32>,
                                                cudaFuncAttributeMaxDynamicSharedMemorySize, int(sm)));
          residues_fused<32><<<grid, 256, sm, stream>>>(d, entries, NB, T, primes, minv, R, plane, e_off, gc);
        }
      RNS_CUDA_CHECK(cudaGetLastError());
    }
  }
  }
}

namespace sdpb_gpu
{
  namespace rns_detail
  {
  namespace
  {
    // Fused CRT (no int32 intermediates): per 16-entry tile and warp,
    //   y_i = q_i (M/m_i)^-1 mod m_i (centred, y = 128 yh + yl) into shared
    //   memory, k = round(sum_i y_i / m_i),
    //   Z (16 x E) = Yh D * 128 + Yl D on INT8 tensor cores (WMMA),
    //   x = sum_j (Z_j - k Mdig_j) 256^j carried in-thread,
    // written as two's complement limbs out[w*out_plane + e0 + e].
    // D: [E][rp] int8 (digit j of M/m_i at j*rp + i), Mdig: [E] int32.
    constexpr int CRT_WARPS = 4;
    __global__ void crt_fused(const uint16_t *__restrict__ Qres, int64_t qplane,
                              int64_t e0, int64_t ne, int r, int rp, int E,
                              const int32_t *__restrict__ primes,
                              const double *__restrict__ minv,
                              const int32_t *__restrict__ inv,
                              const int8_t *__restrict__ D,
                              const int32_t *__restrict__ Mdig,
                              int32_t *__restrict__ out, int64_t out_plane,
                              int *err)
    {
      using namespace nvcuda;
      const int E16 = (E + 15) / 16 * 16, KC = rp / 16;
      extern __shared__ __align__(32) unsigned char smem[];
      int8_t *sD = reinterpret_cast<int8_t *>(smem);              // [KC][E16][16]
      int32_t *sM = reinterpret_cast<int32_t *>(sD + rp * E16);    // [E16]
      unsigned char *wbase = reinterpret_cast<unsigned char *>(sM + E16);
      const size_t per_warp = size_t(2) * 16 * rp + size_t(16) * E16 * 4 + 16 * 16 * 4;
      const int warp = threadIdx.x / 32, lane = threadIdx.x % 32;
      int8_t *sYh = reinterpret_cast<int8_t *>(wbase + warp * per_warp); // [KC][16][16]
      int8_t *sYl = sYh + 16 * rp;
      int32_t *sZ = reinterpret_cast<int32_t *>(sYl + 16 * rp);       // [16][E16]
      int32_t *sT = sZ + 16 * E16;                                     // [16][16]
      for(int q = threadIdx.x; q < rp * E16; q += blockDim.x)
        {
          int kk = q % 16, j = (q / 16) % E16, kc = q / (16 * E16);
          int i = kc * 16 + kk;
          sD[q] = (j < E) ? D[(int64_t)j * rp + i] : 0;
        }
      for(int j = threadIdx.x; j < E16; j += blockDim.x)
        sM[j] = j < E ? Mdig[j] : 0;
      __syncthreads();
      const int64_t tiles = (ne + 15) / 16;
      for(int64_t tile = int64_t(blockIdx.x) * CRT_WARPS + warp; tile < tiles;
          tile += int64_t(gridDim.x) * CRT_WARPS)
        {
          const int64_t t0 = tile * 16;
          for(int q = lane; q < 16 * rp; q += 32)
            {
              int row = q % 16, i = q / 16;
              int32_t yh = 0, yl = 0;
              if(i < r && t0 + row < ne)
                {
                  int32_t m = primes[i];
                  long long qq = Qres[(int64_t)i * qplane + e0 + t0 + row];
                  int32_t y = mod_d(qq * inv[i], m, minv[i]);
                  if(y > m / 2)
                    y -= m;
                  yl = ((y + 64) & 127) - 64;
                  yh = (y - yl) >> 7;
                }
              int off = ((i / 16) * 16 + row) * 16 + i % 16;
              sYh[off] = (int8_t)yh;
              sYl[off] = (int8_t)yl;
            }
          __syncwarp();
          for(int j0 = 0; j0 < E16; j0 += 16)
            {
              wmma::fragment<wmma::accumulator, 16, 16, 16, int> ah, al;
              wmma::fill_fragment(ah, 0);
              wmma::fill_fragment(al, 0);
              for(int kc = 0; kc < KC; ++kc)
                {
                  wmma::fragment<wmma::matrix_a, 16, 16, 16, signed char, wmma::row_major> fh, fl;
                  wmma::fragment<wmma::matrix_b, 16, 16, 16, signed char, wmma::col_major> fb;
                  wmma::load_matrix_sync(fh, sYh + kc * 256, 16);
                  wmma::load_matrix_sync(fl, sYl + kc * 256, 16);
                  wmma::load_matrix_sync(fb, sD + ((size_t)kc * E16 + j0) * 16, 16);
                  wmma::mma_sync(ah, fh, fb, ah);
                  wmma::mma_sync(al, fl, fb, al);
                }
              wmma::store_matrix_sync(sZ + j0, ah, E16, wmma::mem_row_major);
              wmma::store_matrix_sync(sT, al, 16, wmma::mem_row_major);
              __syncwarp();
              for(int q = lane; q < 256; q += 32)
                {
                  int row = q / 16, c = q % 16;
                  sZ[row * E16 + j0 + c] = 128 * sZ[row * E16 + j0 + c] + sT[row * 16 + c];
                }
              __syncwarp();
            }
          if(lane < 16 && t0 + lane < ne)
            {
              const int row = lane;
              double ks = 0;
              for(int i = 0; i < r; ++i)
                {
                  int off = ((i / 16) * 16 + row) * 16 + i % 16;
                  ks += double(128 * int(sYh[off]) + int(sYl[off])) * minv[i];
                }
              const long long k = llrint(ks);
              long long v = 0;
              uint32_t word = 0;
              const int64_t e = e0 + t0 + row;
              for(int j = 0; j < E; ++j)
                {
                  v += (long long)sZ[row * E16 + j] - k * sM[j];
                  uint32_t byte = (uint32_t)(v & 255);
                  v = (v - (long long)byte) >> 8;
                  word |= byte << (8 * (j & 3));
                  if((j & 3) == 3)
                    {
                      out[(int64_t)(j >> 2) * out_plane + e] = (int32_t)word;
                      word = 0;
                    }
                }
              bool neg = int32_t(out[(int64_t)(E / 4 - 1) * out_plane + e]) < 0;
              if(!((v == 0 && !neg) || (v == -1 && neg)))
                atomicExch(err, 1);
            }
          __syncwarp();
        }
    }

    inline void launch_crt_fused(const uint16_t *Qres, int64_t qplane,
                                 int64_t e0, int64_t ne, int r, int rp, int E,
                                 const int32_t *primes, const double *minv,
                                 const int32_t *inv, const int8_t *D,
                                 const int32_t *Mdig, int32_t *out,
                                 int64_t out_plane, int *err,
                                 cudaStream_t stream)
    {
      const int E16 = (E + 15) / 16 * 16;
      const size_t per_warp = size_t(2) * 16 * rp + size_t(16) * E16 * 4 + 16 * 16 * 4;
      const size_t sm = size_t(rp) * E16 + size_t(E16) * 4 + CRT_WARPS * per_warp;
      RNS_CUDA_CHECK(cudaFuncSetAttribute(crt_fused,
                                          cudaFuncAttributeMaxDynamicSharedMemorySize, int(sm)));
      const int64_t tiles = (ne + 15) / 16;
      int grid = int(std::min<int64_t>((tiles + CRT_WARPS - 1) / CRT_WARPS, 132 * 4));
      crt_fused<<<grid, 32 * CRT_WARPS, sm, stream>>>(Qres, qplane, e0, ne, r, rp, E, primes, minv,
                                                     inv, D, Mdig, out, out_plane, err);
      RNS_CUDA_CHECK(cudaGetLastError());
    }
  }
  }
}
