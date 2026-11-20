// GPU backend for SDPB: exact Gram accumulation in residue form.  See
// Gram_Accumulator.hxx for the data layout and the pipeline.
//
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: sdpvibe/rns.py (_residue_kernel: residues as
// an INT8 product of limb bytes with a table of 256^j mod m_i; _crt_dot_kernel /
// crt_carry: CRT as an INT8 product with the digits of M/m_i, then carries),
// sdpvibe/batched.py (_syrk_lower_kernel) and sdpvibe/mpf.py (GramAcc: residues
// accumulated over chunks, one CRT).  Portions adapted from SDPviBe, Copyright
// (c) 2026 Antoine Tilloy, MIT License (docs/GPU_backend_credits.md).
//
// New here: three INT8 products per prime (Karatsuba on r = 128 h + l, with
// 14-bit primes so that h + l fits in int8), cuBLASLt (TMA/wgmma on Hopper)
// for every INT8 product, balanced-digit inputs, the lower triangle cut into
// recursive rectangles.

#include "Gram_Accumulator.hxx"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "rns_kernels.cuh"

#define GPU_CHECK(x) RNS_CUDA_CHECK(x)
#define LT_CHECK(x) RNS_LT_CHECK(x)

namespace sdpb_gpu
{
  namespace
  {
    using namespace rns_detail;

    template <class T> T *dalloc(size_t count, size_t &total)
    {
      T *p = nullptr;
      GPU_CHECK(cudaMalloc(&p, std::max<size_t>(count, 1) * sizeof(T)));
      total += count * sizeof(T);
      return p;
    }

    // Gram block C = R[:, r0:r0+m]^T R[:, c0:c0+n] for the 3g planes (h, l, s
    // of each prime), column-major m x n each.  Q mod m_p += 2^14 HH
    // + 2^7 (SS - HH - LL) + LL, stored in the packed lower triangle.
    __global__ void gram_reduce(const int32_t *__restrict__ C, int m, int n,
                                int g, const int32_t *__restrict__ primes,
                                const double *__restrict__ minv,
                                uint16_t *__restrict__ Qres, int64_t qplane,
                                int N, int r0, int c0)
    {
      int64_t t = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      int64_t mn = (int64_t)m * n;
      if(t >= mn)
        return;
      int a = int(t % m), b = int(t / m);
      int64_t i = r0 + a, j = c0 + b;
      if(i >= N || j >= N || i < j)
        return;
      int64_t pk = j * N - j * (j - 1) / 2 + (i - j);
      for(int p = 0; p < g; ++p)
        {
          long long hh = C[(int64_t)p * mn + t];
          long long ll = C[(int64_t)(g + p) * mn + t];
          long long ss = C[(int64_t)(2 * g + p) * mn + t];
          long long v = 16384LL * hh + 128LL * (ss - hh - ll) + ll;
          int32_t mp = primes[p];
          int32_t r = mod_d(v, mp, minv[p]);
          uint16_t *q = Qres + (int64_t)p * qplane + pk;
          int32_t s = int32_t(*q) + r;
          if(s >= mp)
            s -= mp;
          *q = (uint16_t)s;
        }
    }

  }

  Gram_Accumulator::Gram_Accumulator(int N, int digits, int64_t rows_total,
                                     Gram_Options options, cudaStream_t stream)
      : N_(N), NB_(digits), opt_(options), stream_(stream)
  {
    if(N <= 0 || digits <= 0 || digits % 4 != 0)
      throw std::invalid_argument("Gram_Accumulator: N > 0, digits % 4 == 0");
    Np_ = (N + 31) / 32 * 32;
    g_ = options.prime_group;
    basis_ = Rns_Basis::for_products(digits, rows_total);
    rows_cap_ = 0;
    rows_total_ = rows_total;
    const int r = int(basis_.primes.size()), rp = basis_.num_padded,
              E = basis_.crt_digits;

    LT_CHECK(cublasLtCreate(&lt_));
    GPU_CHECK(cudaMalloc(&workspace_, workspace_bytes_));
    device_bytes_ += workspace_bytes_;

    std::vector<int8_t> Thi, Tlo;
    basis_.residue_tables(digits, Thi, Tlo);
    std::vector<int8_t> T(size_t(rp) * 2 * digits, 0);
    for(int i = 0; i < rp; ++i)
      for(int t = 0; t < digits; ++t)
        {
          T[(size_t(i) * 2 + 0) * digits + t] = Thi[size_t(i) * digits + t];
          T[(size_t(i) * 2 + 1) * digits + t] = Tlo[size_t(i) * digits + t];
        }
    std::vector<int32_t> primes(rp, 1);
    std::vector<double> minv(rp, 1.0);
    for(int i = 0; i < r; ++i)
      {
        primes[i] = basis_.primes[i];
        minv[i] = 1.0 / double(basis_.primes[i]);
      }
    auto inv = basis_.crt_inverses();
    std::vector<int8_t> D;
    std::vector<int32_t> Mdig;
    basis_.crt_tables(D, Mdig);

    d_T_ = dalloc<int8_t>(T.size(), device_bytes_);
    d_primes_ = dalloc<int32_t>(rp, device_bytes_);
    d_minv_ = dalloc<double>(rp, device_bytes_);
    d_inv_ = dalloc<int32_t>(rp, device_bytes_);
    d_D_ = dalloc<int8_t>(D.size(), device_bytes_);
    d_Mdig_ = dalloc<int32_t>(E, device_bytes_);
    GPU_CHECK(cudaMemcpy(d_T_, T.data(), T.size(), cudaMemcpyHostToDevice));
    GPU_CHECK(cudaMemcpy(d_primes_, primes.data(), rp * 4, cudaMemcpyHostToDevice));
    GPU_CHECK(cudaMemcpy(d_minv_, minv.data(), rp * 8, cudaMemcpyHostToDevice));
    GPU_CHECK(cudaMemcpy(d_inv_, inv.data(), rp * 4, cudaMemcpyHostToDevice));
    GPU_CHECK(cudaMemcpy(d_D_, D.data(), D.size(), cudaMemcpyHostToDevice));
    GPU_CHECK(cudaMemcpy(d_Mdig_, Mdig.data(), E * 4, cudaMemcpyHostToDevice));

    // (the int32 intermediate of the residues only without fused residues)
    if(!opt_.fused_residues)
      d_resC_ = dalloc<int32_t>(size_t(2 * g_) * opt_.max_entries_per_residue_batch,
                                device_bytes_);
    d_gram_ = dalloc<int32_t>(size_t(3 * g_) * opt_.max_rect * opt_.max_rect,
                              device_bytes_);
    d_Qres_ = dalloc<uint16_t>(size_t(r) * packed_size(), device_bytes_);
    GPU_CHECK(cudaMemsetAsync(d_Qres_, 0, size_t(r) * packed_size() * 2, stream_));
    const int64_t cb = opt_.crt_batch;
    d_Yh_ = dalloc<int8_t>(size_t(cb) * rp, device_bytes_);
    d_Yl_ = dalloc<int8_t>(size_t(cb) * rp, device_bytes_);
    d_Zh_ = dalloc<int32_t>(size_t(cb) * E, device_bytes_);
    d_K_ = dalloc<int32_t>(cb, device_bytes_);
    d_err_ = dalloc<int>(1, device_bytes_);
    GPU_CHECK(cudaMemsetAsync(d_err_, 0, sizeof(int), stream_));
    GPU_CHECK(cudaEventCreate(&ev0_));
    GPU_CHECK(cudaEventCreate(&ev1_));
  }

  Gram_Accumulator::~Gram_Accumulator()
  {
    for(void *p : {(void *)d_T_, (void *)d_primes_, (void *)d_minv_,
                   (void *)d_inv_, (void *)d_D_, (void *)d_Mdig_,
                   (void *)d_resC_, (void *)d_R_, (void *)d_gram_,
                   (void *)d_Qres_, (void *)d_Yh_, (void *)d_Yl_,
                   (void *)d_Zh_, (void *)d_Zl_, (void *)d_K_, (void *)d_err_,
                   workspace_})
      if(p)
        cudaFree(p);
    cublasLtDestroy(lt_);
    cudaEventDestroy(ev0_);
    cudaEventDestroy(ev1_);
  }

  void Gram_Accumulator::tic()
  {
    if(opt_.timing)
      GPU_CHECK(cudaEventRecord(ev0_, stream_));
  }

  void Gram_Accumulator::toc(const std::string &phase)
  {
    if(!opt_.timing)
      return;
    GPU_CHECK(cudaEventRecord(ev1_, stream_));
    GPU_CHECK(cudaEventSynchronize(ev1_));
    float ms = 0;
    GPU_CHECK(cudaEventElapsedTime(&ms, ev0_, ev1_));
    ms_[phase] += ms;
  }

  void Gram_Accumulator::gemm_tn(int m, int n, int k, const int8_t *A,
                                 int lda, int64_t strideA, const int8_t *B,
                                 int ldb, int64_t strideB, int32_t *C,
                                 int ldc, int64_t strideC, int batch)
  {
    ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_, m, n, k, A,
                       lda, strideA, B, ldb, strideB, C, ldc, strideC, batch);
  }

  void Gram_Accumulator::gram_block(int r0, int r1, int c0, int c1, int i0,
                                    int gc, int64_t rows)
  {
    const int m = r1 - r0, n = c1 - c0;
    const int64_t plane = int64_t(Np_) * rows;
    gemm_tn(m, n, int(rows), d_R_ + int64_t(r0) * rows, int(rows), plane,
            d_R_ + int64_t(c0) * rows, int(rows), plane, d_gram_, m,
            int64_t(m) * n, 3 * gc);
    const int64_t mn = int64_t(m) * n;
    gram_reduce<<<blocks_for(mn), threads, 0, stream_>>>(
      d_gram_, m, n, gc, d_primes_ + i0, d_minv_ + i0,
      d_Qres_ + int64_t(i0) * packed_size(), packed_size(), N_, r0, c0);
    GPU_CHECK(cudaGetLastError());
  }

  // lower triangle of the columns [c0, c1): off-diagonal rectangle, recursion
  void Gram_Accumulator::gram_rect(int c0, int c1, int i0, int gc,
                                   int64_t rows)
  {
    if(c0 >= N_)
      return;
    const int R = opt_.max_rect;
    if(c1 - c0 <= opt_.base)
      {
        gram_block(c0, c1, c0, c1, i0, gc, rows);
        return;
      }
    int half = ((c1 - c0) / 2 + 127) / 128 * 128;
    int h = c0 + half;
    if(h >= c1)
      h = c1 - 32;
    for(int rr = h; rr < c1 && rr < N_; rr += R)
      for(int cc = c0; cc < h; cc += R)
        gram_block(rr, std::min(c1, rr + R), cc, std::min(h, cc + R), i0, gc,
                   rows);
    gram_rect(c0, h, i0, gc, rows);
    gram_rect(h, c1, i0, gc, rows);
  }

  void Gram_Accumulator::add(const int8_t *d_digits, int64_t rows)
  {
    const int64_t cols_per = std::max<int64_t>(
      1, opt_.max_entries_per_residue_batch / std::max<int64_t>(rows, 1));
    add_residues(rows, [&](int i0, int gc, int8_t *R, int64_t plane) {
      if(opt_.fused_residues)
        {
          launch_residues_fused(d_digits, int64_t(N_) * rows, NB_,
                                d_T_ + int64_t(i0) * 2 * NB_, d_primes_ + i0,
                                d_minv_ + i0, R, plane, 0, gc, stream_);
          return;
        }
      for(int64_t j0 = 0; j0 < N_; j0 += cols_per)
        {
          const int64_t nc = std::min<int64_t>(cols_per, N_ - j0);
          const int64_t m_ent = nc * rows;
          gemm_tn(int(m_ent), 2 * gc, NB_, d_digits + j0 * rows * NB_, NB_, 0,
                  d_T_ + int64_t(i0) * 2 * NB_, NB_, 0, d_resC_, int(m_ent),
                  0, 1);
          residue_finish<<<blocks_for(m_ent), threads, 0, stream_>>>(
            d_resC_, m_ent, gc, d_primes_ + i0, d_minv_ + i0, R, plane,
            j0 * rows);
          GPU_CHECK(cudaGetLastError());
        }
    });
  }

  void Gram_Accumulator::add_residues(int64_t rows, const Residue_Fill &fill)
  {
    if(rows <= 0)
      return;
    if(rows % 16 != 0)
      throw std::invalid_argument("Gram_Accumulator::add: rows % 16 != 0 "
                                  "(pad the chunk with zero rows)");
    if(rows > (int64_t(1) << 17))
      throw std::invalid_argument("Gram_Accumulator::add: at most 2^17 rows "
                                  "per chunk (int32 sums of 127^2 products)");
    if(rows_added_ + rows > rows_total_)
      throw std::invalid_argument("Gram_Accumulator::add: more rows than the "
                                  "bound the primes were chosen for");
    const int r = int(basis_.primes.size());
    const int64_t plane = int64_t(Np_) * rows;
    if(rows > rows_cap_)
      {
        if(d_R_)
          {
            cudaFree(d_R_);
            device_bytes_ -= size_t(3 * g_) * Np_ * rows_cap_;
          }
        d_R_ = dalloc<int8_t>(size_t(3 * g_) * Np_ * rows, device_bytes_);
        rows_cap_ = rows;
      }
    for(int i0 = 0; i0 < r; i0 += g_)
      {
        const int gc = std::min(g_, r - i0);
        tic();
        if(Np_ > N_)
          GPU_CHECK(cudaMemset2DAsync(d_R_ + int64_t(N_) * rows, plane, 0,
                                      size_t(Np_ - N_) * rows, 3 * gc,
                                      stream_));
        fill(i0, gc, d_R_, plane);
        toc("residues");
        tic();
        gram_rect(0, Np_, i0, gc, rows);
        toc("gram GEMM + reduce");
      }
    rows_added_ += rows;
  }

  void Gram_Accumulator::finish(int32_t *d_out)
  {
    const int r = int(basis_.primes.size()), rp = basis_.num_padded,
              E = basis_.crt_digits;
    const int64_t total = packed_size(), cb = opt_.crt_batch;
    tic();
    if(opt_.fused_crt)
      launch_crt_fused(d_Qres_, total, 0, total, r, rp, E, d_primes_, d_minv_,
                       d_inv_, d_D_, d_Mdig_, d_out, total, d_err_, stream_);
    for(int64_t e0 = 0; e0 < total && !opt_.fused_crt; e0 += cb)
      {
        const int64_t ne = std::min(cb, total - e0);
        const int64_t ne_pad = (ne + 3) / 4 * 4;
        crt_pre<<<blocks_for(ne_pad), threads, 0, stream_>>>(
          d_Qres_, total, e0, ne, ne_pad, r, rp, d_primes_, d_minv_, d_inv_,
          d_Yh_, d_Yl_, d_K_);
        GPU_CHECK(cudaGetLastError());
        ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_,
                           int(ne_pad), E, rp, d_Yl_, rp, 0, d_D_, rp, 0,
                           d_Zh_, int(ne_pad), 0, 1, 1, 0);
        ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_,
                           int(ne_pad), E, rp, d_Yh_, rp, 0, d_D_, rp, 0,
                           d_Zh_, int(ne_pad), 0, 1, 128, 1);
        crt_carry<<<blocks_for(ne), threads, 0, stream_>>>(
          d_Zh_, nullptr, d_K_, d_Mdig_, ne, ne_pad, E, d_out, total, e0,
          d_err_);
        GPU_CHECK(cudaGetLastError());
      }
    int err = 0;
    GPU_CHECK(cudaMemcpyAsync(&err, d_err_, sizeof(int), cudaMemcpyDeviceToHost,
                              stream_));
    GPU_CHECK(cudaStreamSynchronize(stream_));
    toc("CRT");
    if(err)
      throw std::runtime_error("Gram_Accumulator::finish: CRT carry check "
                               "failed (result out of range)");
    GPU_CHECK(cudaMemsetAsync(d_Qres_, 0, size_t(r) * total * 2, stream_));
    rows_added_ = 0;
  }
}
