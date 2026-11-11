// GPU backend for SDPB: exact batched products in residue arithmetic.  See
// Rns_Gemm.hxx.  Based on SDPviBe by Antoine Tilloy
// (https://github.com/AntoineTilloy/SDPviBe), the prototype of this backend:
// sdpvibe/batched.py (bgemm_nt, _bgemm_kernel) and sdpvibe/rns.py.  Portions
// adapted from SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).

#include "Rns_Gemm.hxx"
#include "rns_kernels.cuh"

#include <algorithm>
#include <cmath>

namespace sdpb_gpu
{
  namespace
  {
    using namespace rns_detail;

    template <class T> T *dalloc(size_t count, size_t &total)
    {
      T *p = nullptr;
      RNS_CUDA_CHECK(cudaMalloc(&p, std::max<size_t>(count, 1) * sizeof(T)));
      total += count * sizeof(T);
      return p;
    }

    // products of 3g planes x batch blocks (prod[(X*g + p)*cplane + e], with
    // cplane = batch*m*pc) -> residues of C mod m_p (Karatsuba)
    __global__ void gemm_reduce(const int32_t *__restrict__ prod,
                                int64_t cplane, int g,
                                const int32_t *__restrict__ primes,
                                const double *__restrict__ minv,
                                uint16_t *__restrict__ Cres)
    {
      int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(e >= cplane)
        return;
      for(int p = 0; p < g; ++p)
        {
          long long hh = prod[(int64_t)p * cplane + e];
          long long ll = prod[(int64_t)(g + p) * cplane + e];
          long long ss = prod[(int64_t)(2 * g + p) * cplane + e];
          long long v = 16384LL * hh + 128LL * (ss - hh - ll) + ll;
          Cres[(int64_t)p * cplane + e] = (uint16_t)mod_d(v, primes[p], minv[p]);
        }
    }

    // gemm_reduce, and the CRT's inputs in the same pass:
    // y_i = c_i (M/m_i)^-1 mod m_i, centred, as int8 halves y = 128 h + l
    // written straight into the CRT rows (l at Y[e*ys + i], h at Y[e*ys + rp +
    // i]), and the quotient estimate Kf[e] = sum_i y_i / m_i accumulated over
    // the prime groups in float32 (|error| < 1e-3 against the margin 1/4 that
    // M > 4|C| leaves).  Replaces the uint16 residues of C and crt_pre.
    __global__ void gemm_reduce_y(const int32_t *__restrict__ prod,
                                  int64_t cplane, int g, int i0,
                                  const int32_t *__restrict__ primes,
                                  const double *__restrict__ minv,
                                  const int32_t *__restrict__ inv,
                                  const float *__restrict__ minvf, int rp,
                                  int8_t *__restrict__ Y, int64_t ys,
                                  float *__restrict__ Kf, bool first)
    {
      // up to 32 primes per group: the entry's bytes are assembled in
      // registers and stored as whole 16-byte pieces (rows start at multiples
      // of 16: ys = 2 rp, rp % 16 == 0, i0 % 32 == 0)
      int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
      if(e >= cplane)
        return;
      float ks = first ? 0.f : Kf[e];
      uint32_t wl[8], wh[8];
#pragma unroll
      for(int w = 0; w < 8; ++w)
        wl[w] = wh[w] = 0;
      const int32_t *ph = prod + e, *pl = ph + g * cplane, *ps = pl + g * cplane;
#pragma unroll
      for(int q = 0; q < 32; ++q)
        if(q < g)
          {
            const int i = i0 + q;
            const int32_t mid = *ps - *ph - *pl;
            const long long v
              = 16384LL * *ph + 128LL * mid + *pl;
            ph += cplane;
            pl += cplane;
            ps += cplane;
            const int32_t m = primes[i];
            const int32_t c = mod_low(v, m, minv[i]);
            // y = c inv mod m, c inv < 2^28: float quotient, exact remainder;
            // the quotient is off by at most one (relative error 2^-23, value
            // below 2^15)
            const uint32_t t = uint32_t(c) * uint32_t(inv[i]);
            int32_t y = int32_t(t) - __float2int_rz(__uint2float_rn(t) * minvf[i]) * m;
            y += (y >> 31) & m;
            y = y >= m ? y - m : y;
            y = y > m / 2 ? y - m : y;
            ks = fmaf(float(y), minvf[i], ks);
            const int32_t lo = ((y + 64) & 127) - 64, hi = (y - lo) >> 7;
            wl[q >> 2] |= uint32_t(uint8_t(int8_t(lo))) << (8 * (q & 3));
            wh[q >> 2] |= uint32_t(uint8_t(int8_t(hi))) << (8 * (q & 3));
          }
      int8_t *yl = Y + e * ys + i0, *yh = Y + e * ys + rp + i0;
      const int nw = (g + 3) / 4;
#pragma unroll
      for(int w = 0; w < 8; w += 4)
        if(w + 4 <= nw)
          {
            *reinterpret_cast<uint4 *>(yl + 4 * w)
              = make_uint4(wl[w], wl[w + 1], wl[w + 2], wl[w + 3]);
            *reinterpret_cast<uint4 *>(yh + 4 * w)
              = make_uint4(wh[w], wh[w + 1], wh[w + 2], wh[w + 3]);
          }
        else
          {
#pragma unroll
            for(int u = 0; u < 4; ++u)
              if(w + u < nw)
                {
                  *reinterpret_cast<uint32_t *>(yl + 4 * (w + u)) = wl[w + u];
                  *reinterpret_cast<uint32_t *>(yh + 4 * (w + u)) = wh[w + u];
                }
          }
      Kf[e] = ks;
    }

    std::vector<int8_t> interleave(const std::vector<int8_t> &hi,
                                   const std::vector<int8_t> &lo, int rp,
                                   int D)
    {
      std::vector<int8_t> T(size_t(rp) * 2 * D, 0);
      for(int i = 0; i < rp; ++i)
        for(int t = 0; t < D; ++t)
          {
            T[(size_t(i) * 2 + 0) * D + t] = hi[size_t(i) * D + t];
            T[(size_t(i) * 2 + 1) * D + t] = lo[size_t(i) * D + t];
          }
      return T;
    }
  }

  Rns_Gemm::Rns_Gemm(int m, int n, int k, int batch, int digitsA,
                     int digitsB, Gram_Options options, int panel,
                     cudaStream_t stream)
      : m_(m), n_(n), k_(k), batch_(batch), DA_(digitsA), DB_(digitsB),
        g_(options.prime_group), np_(std::min(panel, n)), opt_(options),
        stream_(stream)
  {
    if(m % 4 || k % 16 || digitsA % 4 || digitsB % 4)
      throw std::invalid_argument("Rns_Gemm: m % 4, k % 16, digits % 4");
    // |C| < k 2^(8(DA+DB) - 2); M > 4 |C|
    basis_ = Rns_Basis::for_bound(
      8 * (digitsA + digitsB) + int(std::ceil(std::log2(double(k)))));
    const int r = int(basis_.primes.size()), rp = basis_.num_padded,
              E = basis_.crt_digits;
    // a group is at most all the primes: the buffers are sized by it
    g_ = std::min(g_, r);
    RNS_LT_CHECK(cublasLtCreate(&lt_));
    // cuBLASLt's scratch: one per product unless shared (64 MiB each, and
    // adaptive digits make many products of one shape)
    if(!opt_.workspace)
      {
        RNS_CUDA_CHECK(cudaMalloc(&workspace_, workspace_bytes_));
        device_bytes_ += workspace_bytes_;
      }
    std::vector<int8_t> hi, lo;
    basis_.residue_tables(digitsA, hi, lo);
    auto TA = interleave(hi, lo, rp, digitsA);
    basis_.residue_tables(digitsB, hi, lo);
    auto TB = interleave(hi, lo, rp, digitsB);
    std::vector<int32_t> primes(rp, 1);
    std::vector<double> minv(rp, 1.0);
    for(int i = 0; i < r; ++i)
      {
        primes[i] = basis_.primes[i];
        minv[i] = 1.0 / primes[i];
      }
    auto inv = basis_.crt_inverses();
    std::vector<int8_t> D;
    std::vector<int32_t> Mdig;
    basis_.crt_tables(D, Mdig);
    d_TA_ = dalloc<int8_t>(TA.size(), device_bytes_);
    d_TB_ = dalloc<int8_t>(TB.size(), device_bytes_);
    d_primes_ = dalloc<int32_t>(rp, device_bytes_);
    d_minv_ = dalloc<double>(rp, device_bytes_);
    d_inv_ = dalloc<int32_t>(rp, device_bytes_);
    d_D_ = dalloc<int8_t>(D.size(), device_bytes_);
    d_Mdig_ = dalloc<int32_t>(E, device_bytes_);
    RNS_CUDA_CHECK(cudaMemcpy(d_TA_, TA.data(), TA.size(), cudaMemcpyHostToDevice));
    RNS_CUDA_CHECK(cudaMemcpy(d_TB_, TB.data(), TB.size(), cudaMemcpyHostToDevice));
    RNS_CUDA_CHECK(cudaMemcpy(d_primes_, primes.data(), rp * 4, cudaMemcpyHostToDevice));
    RNS_CUDA_CHECK(cudaMemcpy(d_minv_, minv.data(), rp * 8, cudaMemcpyHostToDevice));
    RNS_CUDA_CHECK(cudaMemcpy(d_inv_, inv.data(), rp * 4, cudaMemcpyHostToDevice));
    RNS_CUDA_CHECK(cudaMemcpy(d_D_, D.data(), D.size(), cudaMemcpyHostToDevice));
    RNS_CUDA_CHECK(cudaMemcpy(d_Mdig_, Mdig.data(), E * 4, cudaMemcpyHostToDevice));

    const int64_t cplane = int64_t(batch) * m * np_;
    const bool shared = opt_.workspace != nullptr;
    // (the int32 intermediate of residues() is used only without fused residues
    // when the digits come from the caller; the panel's own copy of A's
    // residues only without the cache)
    if(!opt_.fused_residues)
      d_resC_ = dalloc<int32_t>(size_t(2 * g_) * opt_.max_entries_per_residue_batch,
                                device_bytes_);
    if(!shared && !opt_.cache_a)
      d_RA_ = dalloc<int8_t>(size_t(3 * g_) * batch * k * m, device_bytes_);
    if(!shared)
      {
        d_RB_ = dalloc<int8_t>(size_t(3 * g_) * batch * k * np_, device_bytes_);
        d_prod_ = dalloc<int32_t>(size_t(3 * g_) * cplane, device_bytes_);
      }
    reduce_crt_ = opt_.reduce_to_crt && !opt_.fused_crt && g_ <= 32;
    if(reduce_crt_)
      {
        // CRT rows for the whole panel (+ 4 rows: the products read entries
        // padded to a multiple of 4). Rows' padding primes are multiplied by the
        // zero rows of the tables, so whatever those bytes hold does not count.
        if(!shared)
          {
            d_Yall_ = dalloc<int8_t>(size_t(cplane + 4) * 2 * rp, device_bytes_);
            d_Kf_ = dalloc<float>(cplane, device_bytes_);
          }
        std::vector<float> mf(rp, 1.0f);
        for(int i = 0; i < r; ++i)
          mf[i] = float(1.0 / primes[i]);
        d_minvf_ = dalloc<float>(rp, device_bytes_);
        RNS_CUDA_CHECK(cudaMemcpy(d_minvf_, mf.data(), rp * 4, cudaMemcpyHostToDevice));
      }
    else if(!shared)
      d_Cres_ = dalloc<uint16_t>(size_t(r) * cplane, device_bytes_);
    // CRT buffers: a batch is at most one panel's entries (padded to 4)
    const int64_t cb = std::min<int64_t>(opt_.crt_batch, (cplane + 3) / 4 * 4);
    if(opt_.crt_base128 && !opt_.fused_crt)
      {
        std::vector<int8_t> D2;
        std::vector<int32_t> Mdig2;
        basis_.crt_tables_128(D2, Mdig2);
        E2_ = basis_.crt_digits_128();
        if(opt_.crt_guard_bits != 0)
          {
            // the digits below j0 sum to at most 264 * 64 * 64 * 128^j0 / 127
            // (each of the 2 r products |l|, |h| <= 64 by a balanced digit
            // <= 64; the k M term is smaller): below 2^13.1 128^j0, so
            // 7 j0 + 14 <= 8 D - 2 - guard
            const int D = std::min(digitsA, digitsB);
            j0_ = std::max(0, (8 * D - 2 - opt_.crt_guard_bits - 14) / 7);
            Wt_ = (7 * (E2_ - j0_) + 31) / 32 + 1;
          }
        d_D2_ = dalloc<int8_t>(D2.size(), device_bytes_);
        d_Mdig2_ = dalloc<int32_t>(E2_, device_bytes_);
        RNS_CUDA_CHECK(cudaMemcpy(d_D2_, D2.data(), D2.size(), cudaMemcpyHostToDevice));
        RNS_CUDA_CHECK(cudaMemcpy(d_Mdig2_, Mdig2.data(), E2_ * 4, cudaMemcpyHostToDevice));
        if(!shared)
          {
            if(!reduce_crt_)
              d_Y2_ = dalloc<int8_t>(size_t(cb) * 2 * rp, device_bytes_);
            d_Z2_ = dalloc<int32_t>(size_t(cb) * (E2_ - j0_), device_bytes_);
          }
      }
    else
      {
        d_Yh_ = dalloc<int8_t>(size_t(cb) * rp, device_bytes_);
        d_Yl_ = dalloc<int8_t>(size_t(cb) * rp, device_bytes_);
        d_Zh_ = dalloc<int32_t>(size_t(cb) * E, device_bytes_);
      }
    if(!shared)
      d_K_ = dalloc<int32_t>(cb, device_bytes_);
    d_err_ = dalloc<int>(1, device_bytes_);
    RNS_CUDA_CHECK(cudaMemsetAsync(d_err_, 0, sizeof(int), stream_));
    RNS_CUDA_CHECK(cudaEventCreate(&ev0_));
    RNS_CUDA_CHECK(cudaEventCreate(&ev1_));
  }

  Rns_Workspace::~Rns_Workspace()
  {
    for(void *p : p_)
      if(p)
        cudaFree(p);
  }

  void *Rns_Workspace::get(Buffer b, size_t bytes)
  {
    if(n_[b] < bytes)
      {
        RNS_CUDA_CHECK(cudaStreamSynchronize(stream_));
        if(p_[b])
          cudaFree(p_[b]);
        p_[b] = nullptr;
        n_[b] = 0;
        RNS_CUDA_CHECK(cudaMalloc(&p_[b], bytes));
        n_[b] = bytes;
      }
    return p_[b];
  }

  size_t Rns_Workspace::bytes() const
  {
    size_t total = 0;
    for(size_t n : n_)
      total += n;
    return total;
  }

  void Rns_Workspace::release()
  {
    RNS_CUDA_CHECK(cudaStreamSynchronize(stream_));
    for(int b = 0; b < count; ++b)
      {
        if(p_[b])
          cudaFree(p_[b]);
        p_[b] = nullptr;
        n_[b] = 0;
      }
  }

  void Rns_Gemm::bind_workspace()
  {
    Rns_Workspace *w = opt_.workspace;
    if(!w)
      return;
    const int r = int(basis_.primes.size()), rp = basis_.num_padded;
    const int64_t cplane = int64_t(batch_) * m_ * np_;
    const int64_t cb = std::min<int64_t>(opt_.crt_batch, (cplane + 3) / 4 * 4);
    const size_t planeA = size_t(batch_) * k_ * m_;
    using W = Rns_Workspace;
    workspace_ = w->get(W::lt, workspace_bytes_);
    if(opt_.cache_a)
      d_RAall_ = static_cast<int8_t *>(w->get(W::RA_all, size_t(3) * r * planeA));
    else
      d_RA_ = static_cast<int8_t *>(w->get(W::RA, size_t(3 * g_) * planeA));
    d_RB_ = static_cast<int8_t *>(
      w->get(W::RB, size_t(3 * g_) * batch_ * k_ * np_));
    d_prod_ = static_cast<int32_t *>(w->get(W::prod, size_t(3 * g_) * cplane * 4));
    if(reduce_crt_)
      {
        d_Yall_ = static_cast<int8_t *>(
          w->get(W::Y_all, size_t(cplane + 4) * 2 * rp));
        d_Kf_ = static_cast<float *>(w->get(W::Kf, size_t(cplane) * 4));
      }
    else
      {
        d_Cres_ = static_cast<uint16_t *>(w->get(W::Cres, size_t(r) * cplane * 2));
        if(E2_)
          d_Y2_ = static_cast<int8_t *>(w->get(W::Y2, size_t(cb) * 2 * rp));
      }
    if(E2_)
      d_Z2_ = static_cast<int32_t *>(
        w->get(W::Z2, size_t(cb) * (E2_ - j0_) * 4));
    d_K_ = static_cast<int32_t *>(w->get(W::K, size_t(cb) * 4));
  }

  Rns_Gemm::~Rns_Gemm()
  {
    if(opt_.workspace)
      d_RA_ = d_RB_ = d_RAall_ = d_Yall_ = d_Y2_ = nullptr,
      d_prod_ = d_Z2_ = d_K_ = nullptr, d_Kf_ = nullptr, d_Cres_ = nullptr,
      workspace_ = nullptr;
    for(void *p :
        {(void *)d_TA_, (void *)d_TB_, (void *)d_primes_, (void *)d_minv_,
         (void *)d_inv_, (void *)d_D_, (void *)d_Mdig_, (void *)d_resC_,
         (void *)d_RA_, (void *)d_RB_, (void *)d_prod_, (void *)d_Cres_,
         (void *)d_Yh_, (void *)d_Yl_, (void *)d_Zh_, (void *)d_Zl_,
         (void *)d_K_, (void *)d_err_, (void *)d_RAall_, (void *)d_D2_,
         (void *)d_Y2_, (void *)d_Mdig2_, (void *)d_Z2_, (void *)d_Yall_,
         (void *)d_Kf_, (void *)d_minvf_, workspace_})
      if(p)
        cudaFree(p);
    cublasLtDestroy(lt_);
    cudaEventDestroy(ev0_);
    cudaEventDestroy(ev1_);
  }

  void Rns_Gemm::tic()
  {
    if(opt_.timing)
      RNS_CUDA_CHECK(cudaEventRecord(ev0_, stream_));
  }

  void Rns_Gemm::toc(const std::string &phase)
  {
    if(!opt_.timing)
      return;
    RNS_CUDA_CHECK(cudaEventRecord(ev1_, stream_));
    RNS_CUDA_CHECK(cudaEventSynchronize(ev1_));
    float ms = 0;
    RNS_CUDA_CHECK(cudaEventElapsedTime(&ms, ev0_, ev1_));
    ms_[phase] += ms;
  }

  // residues of batch_ blocks of ent_per_block contiguous entries each (block
  // b at d + b*stride*digits) into planes R[(X*gc + p)*plane + b*ent + e]
  void Rns_Gemm::residues(const int8_t *d, int64_t ent, int64_t stride,
                          int digits, const int8_t *T, int i0, int gc,
                          int8_t *R, int64_t plane)
  {
    const int64_t maxe = opt_.max_entries_per_residue_batch;
    if(opt_.fused_residues)
      {
        for(int b = 0; b < batch_; ++b)
          launch_residues_fused(d + b * stride * digits, ent, digits,
                                T + int64_t(i0) * 2 * digits, d_primes_ + i0,
                                d_minv_ + i0, R, plane, b * ent, gc, stream_);
        return;
      }
    for(int b = 0; b < batch_; ++b)
      for(int64_t e0 = 0; e0 < ent; e0 += maxe)
        {
          const int64_t ne = std::min(maxe, ent - e0);
          ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_,
                             int(ne), 2 * gc, digits,
                             d + (b * stride + e0) * digits, digits, 0,
                             T + int64_t(i0) * 2 * digits, digits, 0, d_resC_,
                             int(ne), 0, 1);
          residue_finish<<<blocks_for(ne), threads, 0, stream_>>>(
            d_resC_, ne, gc, d_primes_ + i0, d_minv_ + i0, R, plane,
            b * ent + e0);
          RNS_CUDA_CHECK(cudaGetLastError());
        }
  }

  void Rns_Gemm::multiply(const int8_t *dA, const int8_t *dB, int32_t *d_out,
                          int64_t strideA, int64_t strideB)
  {
    if(strideA < 0)
      strideA = int64_t(k_) * m_;
    if(strideB < 0)
      strideB = int64_t(k_) * n_;
    multiply_with(
      [&](int i0, int gc, int8_t *RA, int64_t planeA) {
        residues(dA, int64_t(k_) * m_, strideA, DA_, d_TA_, i0, gc, RA, planeA);
      },
      [&](int i0, int gc, int j0, int pc, int8_t *RB, int64_t planeB) {
        residues(dB + int64_t(j0) * k_ * DB_, int64_t(k_) * pc, strideB, DB_,
                 d_TB_, i0, gc, RB, planeB);
      },
      d_out);
  }

  void Rns_Gemm::multiply_with(const Fill_A &fillA, const Fill_B &fillB,
                               int32_t *d_out, const Consume &consume,
                               int triangles)
  {
    const bool same = (triangles & Product_Same_Operands) && opt_.cache_a
                      && m_ == n_ && DA_ == DB_;
    const bool lower = triangles & Product_Lower_Only;
    bind_workspace();
    triangles &= ~(Product_Same_Operands | Product_Lower_Only | Product_Negate);
    if((triangles != 0 && triangles != 1 && triangles != 2 && triangles != 4)
       || (triangles == 1 && m_ != k_)
       || (triangles >= 2 && n_ != k_) || (lower && (m_ != n_ || d_out)))
      throw std::invalid_argument("Rns_Gemm: triangular shape or flag");
    const int r = int(basis_.primes.size());
    const int64_t planeA = int64_t(batch_) * k_ * m_;
    // lower triangle only: narrower panels leave out more of the upper triangle
    // (four panels omit 3/8 of C; each costs one more launch per prime group)
    const int width
      = lower ? std::min(np_, std::max(128, ((m_ + 3) / 4 + 15) / 16 * 16)) : np_;
    for(int j0 = 0; j0 < n_; j0 += width)
      {
        const int pc = std::min(width, n_ - j0);
        // rows of C wanted in this panel
        const int row0 = lower ? j0 : 0, rows = m_ - row0;
        const int64_t cplane = int64_t(batch_) * rows * pc;
        const int64_t planeB = int64_t(batch_) * k_ * pc;
        for(int i0 = 0; i0 < r; i0 += g_)
          {
            const int gc = std::min(g_, r - i0);
            tic();
            // A's residues: once for all panels when cached (group i0's 3 gc
            // planes start at plane 3 i0 of the cache), else per panel
            int8_t *RA = d_RA_;
            if(borrowed_RA_)
              RA = const_cast<int8_t *>(borrowed_RA_) + int64_t(3) * i0 * planeA;
            else if(opt_.cache_a)
              {
                if(!d_RAall_)
                  d_RAall_ = dalloc<int8_t>(size_t(3) * r * planeA, device_bytes_);
                RA = d_RAall_ + int64_t(3) * i0 * planeA;
                if(j0 == 0)
                  fillA(i0, gc, RA, planeA);
              }
            else
              fillA(i0, gc, RA, planeA);
            toc("residues A");
            tic();
            // B's residues, or the columns of the same matrix in A's cache
            int8_t *RB = same ? RA + int64_t(j0) * k_ : d_RB_;
            const int64_t strideB = same ? int64_t(k_) * m_ : int64_t(k_) * pc;
            if(!same)
              fillB(i0, gc, j0, pc, d_RB_, planeB);
            toc("residues B");
            tic();
            // with a triangular operand, known zeros are not multiplied:
            // tiles of rows if it is A, the inner range of k if it is B
            const int row_tile = triangles == 1 ? np_ : rows;
            for(int row = row0; row < m_; row += row_tile)
              {
                const int nr = std::min(row_tile, m_ - row);
                const int begin = triangles == 2 ? j0 : 0;
                const int end = triangles == 1   ? row + nr
                                : triangles == 4 ? j0 + pc
                                                 : k_;
                ops_ += lt_gemm_tn(
                  lt_, workspace_, workspace_bytes_, stream_, nr, pc, end - begin,
                  RA + int64_t(row) * k_ + begin, k_, int64_t(k_) * m_,
                  RB + begin, k_, strideB, d_prod_ + (row - row0), rows,
                  int64_t(rows) * pc, 3 * gc * batch_);
              }
            toc("block GEMMs");
            reduce_group(cplane, i0, gc);
          }
        crt_panel(cplane, j0, pc, d_out, consume);
      }
    borrowed_RA_ = nullptr;
    if(!opt_.defer_crt_check)
      check();
  }

  void Rns_Gemm::borrow_a_residues(const Rns_Gemm &other)
  {
    if(other.d_RAall_ && other.m_ == m_ && other.k_ == k_ && other.batch_ == batch_
       && other.DA_ == DA_ && other.g_ == g_
       && other.basis_.primes.size() == basis_.primes.size())
      borrowed_RA_ = other.d_RAall_;
  }

  void Rns_Gemm::reduce_group(int64_t cplane, int i0, int gc)
  {
    const int rp = basis_.num_padded;
    tic();
    if(d_Yall_)
      gemm_reduce_y<<<blocks_for(cplane), threads, 0, stream_>>>(
        d_prod_, cplane, gc, i0, d_primes_, d_minv_, d_inv_, d_minvf_,
        rp, d_Yall_, 2 * rp, d_Kf_, i0 == 0);
    else
      gemm_reduce<<<blocks_for(cplane), threads, 0, stream_>>>(
        d_prod_, cplane, gc, d_primes_ + i0, d_minv_ + i0,
        d_Cres_ + int64_t(i0) * cplane);
    RNS_CUDA_CHECK(cudaGetLastError());
    toc("reduce");
  }

  void Rns_Gemm::crt_panel(int64_t cplane, int j0, int pc, int32_t *d_out,
                           const Consume &consume)
  {
    const int r = int(basis_.primes.size()), rp = basis_.num_padded,
              E = basis_.crt_digits, W = limbs_out();
    // CRT of the panel, written at out[w*plane + b*m*n + (j0 + jj)*m + i]
    tic();
    const int64_t cb = opt_.crt_batch;
    const bool direct = (pc == n_) && d_out;
    int32_t *tmp = d_out;
    const bool shared_tmp = !direct && opt_.workspace;
    if(shared_tmp)
      tmp = static_cast<int32_t *>(opt_.workspace->get(
        Rns_Workspace::exact, size_t(W) * int64_t(batch_) * m_ * np_ * 4));
    else if(!direct)
      RNS_CUDA_CHECK(malloc_async_retry(reinterpret_cast<void **>(&tmp),
                                        size_t(W) * cplane * 4, stream_));
    if(opt_.fused_crt)
      launch_crt_fused(d_Cres_, cplane, 0, cplane, r, rp, E, d_primes_,
                       d_minv_, d_inv_, d_D_, d_Mdig_, tmp, cplane, d_err_,
                       stream_);
    for(int64_t e0 = 0; e0 < cplane && E2_; e0 += cb)
      {
        // one product: Z = [l | h] D2 (base-128 digits; Rns_Basis::crt_tables_128)
        const int64_t ne = std::min(cb, cplane - e0);
        const int64_t ne_pad = (ne + 3) / 4 * 4;
        const int8_t *Y = d_Yall_ ? d_Yall_ + e0 * 2 * rp : d_Y2_;
        if(!d_Yall_)
          {
            crt_pre<<<blocks_for(ne_pad), threads, 0, stream_>>>(
              d_Cres_, cplane, e0, ne, ne_pad, r, rp, d_primes_, d_minv_,
              d_inv_, d_Y2_ + rp, d_Y2_, d_K_, 2 * rp);
            RNS_CUDA_CHECK(cudaGetLastError());
          }
        toc("CRT 1 digits of y");
        tic();
        ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_,
                           int(ne_pad), E2_ - j0_, 2 * rp, Y, 2 * rp, 0,
                           d_D2_ + int64_t(j0_) * 2 * rp, 2 * rp, 0, d_Z2_,
                           int(ne_pad), 0, 1, 1, 0);
        toc("CRT 2 products");
        tic();
        crt_carry_128<<<blocks_for(ne), threads, 0, stream_>>>(
          d_Z2_, d_K_, d_Mdig2_ + j0_, ne, ne_pad, E2_ - j0_, W, tmp, cplane, e0,
          d_err_, d_Yall_ ? d_Kf_ + e0 : nullptr);
        RNS_CUDA_CHECK(cudaGetLastError());
        toc("CRT 3 carries");
        tic();
      }
    for(int64_t e0 = 0; e0 < cplane && !opt_.fused_crt && !E2_; e0 += cb)
      {
        const int64_t ne = std::min(cb, cplane - e0);
        const int64_t ne_pad = (ne + 3) / 4 * 4;
        const int8_t *Yl = d_Yall_ ? d_Yall_ + e0 * 2 * rp : d_Yl_,
                     *Yh = d_Yall_ ? d_Yall_ + e0 * 2 * rp + rp : d_Yh_;
        const int ldy = d_Yall_ ? 2 * rp : rp;
        if(!d_Yall_)
          {
            crt_pre<<<blocks_for(ne_pad), threads, 0, stream_>>>(
              d_Cres_, cplane, e0, ne, ne_pad, r, rp, d_primes_, d_minv_,
              d_inv_, d_Yh_, d_Yl_, d_K_);
            RNS_CUDA_CHECK(cudaGetLastError());
          }
        toc("CRT 1 digits of y");
        tic();
        // Z = 128 (Yh D) + Yl D accumulated in one int32 buffer
        ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_,
                           int(ne_pad), E, rp, Yl, ldy, 0, d_D_, rp, 0,
                           d_Zh_, int(ne_pad), 0, 1, 1, 0);
        ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_,
                           int(ne_pad), E, rp, Yh, ldy, 0, d_D_, rp, 0,
                           d_Zh_, int(ne_pad), 0, 1, 128, 1);
        toc("CRT 2 products");
        tic();
        crt_carry<<<blocks_for(ne), threads, 0, stream_>>>(
          d_Zh_, nullptr, d_K_, d_Mdig_, ne, ne_pad, E, tmp, cplane, e0,
          d_err_, d_Yall_ ? d_Kf_ + e0 : nullptr);
        RNS_CUDA_CHECK(cudaGetLastError());
        toc("CRT 3 carries");
        tic();
      }
    toc("CRT 0 other");
    tic();
    if(!d_out)
      {
        consume(tmp, cplane, j0, pc);
        if(!shared_tmp)
          RNS_CUDA_CHECK(cudaFreeAsync(tmp, stream_));
      }
    else if(!direct)
      {
        RNS_CUDA_CHECK(cudaMemcpy2DAsync(
          d_out + int64_t(j0) * m_, size_t(m_) * n_ * 4, tmp,
          size_t(m_) * pc * 4, size_t(m_) * pc * 4, size_t(W) * batch_,
          cudaMemcpyDeviceToDevice, stream_));
        if(!shared_tmp)
          RNS_CUDA_CHECK(cudaFreeAsync(tmp, stream_));
      }
    toc("CRT 4 output (rounding)");
  }

  void Rns_Gemm::multiply_panels(const Fill_A &fillA, const Fill_Panel &fillB,
                                 const Consume_Panel &consume,
                                 const std::vector<Product_Panel> &panels)
  {
    const int r = int(basis_.primes.size());
    const int64_t planeA = int64_t(batch_) * k_ * m_;
    bind_workspace();
    bool first = true;
    for(const auto &p : panels)
      {
        const int rows = m_ - p.row0;
        const int64_t cplane = int64_t(batch_) * rows * p.pc;
        const int64_t planeB = int64_t(batch_) * k_ * p.pc;
        for(int i0 = 0; i0 < r; i0 += g_)
          {
            const int gc = std::min(g_, r - i0);
            tic();
            int8_t *RA = d_RA_;
            if(opt_.cache_a)
              {
                if(!d_RAall_)
                  d_RAall_ = dalloc<int8_t>(size_t(3) * r * planeA, device_bytes_);
                RA = d_RAall_ + int64_t(3) * i0 * planeA;
                if(first)
                  fillA(i0, gc, RA, planeA);
              }
            else
              fillA(i0, gc, RA, planeA);
            toc("residues A");
            tic();
            fillB(i0, gc, p, d_RB_, planeB);
            toc("residues B");
            tic();
            ops_ += lt_gemm_tn(lt_, workspace_, workspace_bytes_, stream_, rows,
                               p.pc, p.last - p.first,
                               RA + int64_t(p.row0) * k_ + p.first, k_,
                               int64_t(k_) * m_, d_RB_ + p.first, k_,
                               int64_t(k_) * p.pc, d_prod_, rows,
                               int64_t(rows) * p.pc, 3 * gc * batch_);
            toc("block GEMMs");
            reduce_group(cplane, i0, gc);
          }
        first = false;
        crt_panel(cplane, p.j0, p.pc, nullptr,
                  [&](const int32_t *panel, int64_t c, int, int) {
                    consume(panel, c, p);
                  });
      }
    if(!opt_.defer_crt_check)
      check();
  }

  void Rns_Gemm::check()
  {
    int err = 0;
    RNS_CUDA_CHECK(cudaMemcpyAsync(&err, d_err_, sizeof(int),
                                   cudaMemcpyDeviceToHost, stream_));
    RNS_CUDA_CHECK(cudaStreamSynchronize(stream_));
    if(err)
      throw std::runtime_error("Rns_Gemm: CRT carry check failed");
  }
}
