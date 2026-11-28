// GPU backend for SDPB: the element kernels, matrix-vector products and
// reductions of solver/Mpf_Ops against exact arithmetic in GMP, which shares
// nothing with them.  Exit 0 if every bound holds.
//   mpf_ops_check [--perturb-bits b]
// --perturb-bits b: adds 2^-b relative to the matrix read by GMP, the negative
// control: the gemv checks must then fail for b well above the bound.
#include "../solver/Mpf_Ops.hxx"

#include <gmp.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sdpb_gpu;

namespace
{
  int L = 28, NB = 32 * L - 2; // --L
  constexpr int64_t EZ = MPF_EZ;
  std::mt19937_64 rng(332);
  int failures = 0;

  struct Host_Mat
  {
    int batch;
    int64_t rows, cols;
    std::vector<int32_t> m;
    std::vector<int64_t> e;
    int64_t entries() const { return int64_t(batch) * rows * cols; }
    Host_Mat(int b, int64_t r, int64_t c)
        : batch(b), rows(r), cols(c), m(size_t(L) * b * r * c), e(b * r * c)
    {}
  };

  void set_random(Host_Mat &A, int64_t i, int spread, double zero_prob)
  {
    const int64_t tot = A.entries();
    if(std::uniform_real_distribution<>(0, 1)(rng) < zero_prob)
      {
        for(int t = 0; t < L; ++t)
          A.m[t * tot + i] = 0;
        A.e[i] = EZ;
        return;
      }
    uint32_t w[L];
    for(int t = 0; t < L; ++t)
      w[t] = uint32_t(rng());
    w[L - 1] = (w[L - 1] & 0x1fffffffu) | 0x20000000u; // exactly NB bits
    const bool neg = rng() & 1;
    if(neg)
      {
        unsigned long long c = 1;
        for(int t = 0; t < L; ++t)
          {
            const unsigned long long v = (unsigned long long)(~w[t]) + c;
            w[t] = uint32_t(v);
            c = v >> 32;
          }
      }
    for(int t = 0; t < L; ++t)
      A.m[t * tot + i] = int32_t(w[t]);
    A.e[i] = -NB + int64_t(rng() % (2 * spread + 1)) - spread;
  }

  Host_Mat random_mat(int b, int64_t r, int64_t c, int spread, double zero_prob)
  {
    Host_Mat A(b, r, c);
    for(int64_t i = 0; i < A.entries(); ++i)
      set_random(A, i, spread, zero_prob);
    return A;
  }

  // exact value of entry i (mpf with enough bits for every sum below)
  void value(const int32_t *m, int64_t plane, int64_t i, int64_t e, mpf_t out,
             int limbs = L)
  {
    if(e == EZ)
      {
        mpf_set_ui(out, 0);
        return;
      }
    std::vector<uint32_t> a(limbs);
    for(int t = 0; t < limbs; ++t)
      a[t] = uint32_t(m[t * plane + i]);
    const bool neg = int32_t(a[limbs - 1]) < 0;
    if(neg)
      {
        unsigned long long c = 1;
        for(auto &w : a)
          {
            const unsigned long long v = (unsigned long long)(~w) + c;
            w = uint32_t(v);
            c = v >> 32;
          }
      }
    mpz_t z;
    mpz_init(z);
    mpz_import(z, limbs, -1, 4, 0, 0, a.data());
    mpf_set_z(out, z);
    mpz_clear(z);
    if(e >= 0)
      mpf_mul_2exp(out, out, e);
    else
      mpf_div_2exp(out, out, -e);
    if(neg)
      mpf_neg(out, out);
  }

  double log2_abs(const mpf_t x)
  {
    if(mpf_sgn(x) == 0)
      return -1e30;
    long ex;
    const double d = mpf_get_d_2exp(&ex, x);
    return std::log2(std::fabs(d)) + double(ex);
  }

  struct Dev
  {
    Mpf_Batch b;
    explicit Dev(const Host_Mat &H)
    {
      b.L = L;
      b.batch = H.batch;
      b.rows = H.rows;
      b.cols = H.cols;
      if(cudaMalloc(&b.m, H.m.size() * 4) != cudaSuccess
         || cudaMalloc(&b.e, H.e.size() * 8) != cudaSuccess)
        throw std::runtime_error("cudaMalloc");
      mpf_upload(H.m.data(), H.e.data(), b, 0);
    }
    Dev(int batch, int64_t rows, int64_t cols)
    {
      b.L = L;
      b.batch = batch;
      b.rows = rows;
      b.cols = cols;
      if(cudaMalloc(&b.m, size_t(L) * b.entries() * 4) != cudaSuccess
         || cudaMalloc(&b.e, size_t(b.entries()) * 8) != cudaSuccess)
        throw std::runtime_error("cudaMalloc");
    }
    ~Dev()
    {
      cudaFree(b.m);
      cudaFree(b.e);
    }
    Host_Mat get() const
    {
      Host_Mat H(b.batch, b.rows, b.cols);
      mpf_download(b, H.m.data(), H.e.data(), 0);
      return H;
    }
  };

  void report(const char *name, double worst, double bound)
  {
    const bool ok = worst <= bound;
    printf("%-34s worst 2^%.1f, bound 2^%.1f  %s\n", name, worst, bound,
           ok ? "ok" : "FAIL");
    if(!ok)
      ++failures;
  }

  struct Mpf_Vec
  {
    mpf_t a, b, c, d;
    Mpf_Vec()
    {
      mpf_inits(a, b, c, d, nullptr);
    }
    ~Mpf_Vec() { mpf_clears(a, b, c, d, nullptr); }
  };

  // checks a gemv: index maps as in Gemv_Spec, bound relative to the largest term
  void check_gemv(const char *name, const Host_Mat &A, const Host_Mat &v,
                  const Gemv_Spec &g, double perturb_bits)
  {
    Dev dA(A), dv(v);
    const int64_t oent = int64_t(g.batch) * g.no;
    Dev dO(g.batch, g.no, 1);
    // spec's outputs are b ob + o oo with ob = no here
    mpf_gemv(dA.b, dv.b, dO.b, g, 0);
    cudaStreamSynchronize(0);
    Host_Mat O = dO.get();
    double worst = -1e30;
    Mpf_Vec t;
    mpf_t sum, big;
    mpf_inits(sum, big, nullptr);
    for(int b = 0; b < g.batch; ++b)
      for(int64_t o = 0; o < g.no; ++o)
        {
          mpf_set_ui(sum, 0);
          mpf_set_ui(big, 0);
          for(int64_t tb = 0; tb < g.ntb; ++tb)
            for(int64_t ti = 0; ti < g.nti; ++ti)
              {
                const int64_t ia = b * g.ab + o * g.ao + tb * g.atb + ti * g.ati;
                const int64_t iv = b * g.vb + o * g.vo + tb * g.vtb + ti * g.vti;
                value(A.m.data(), A.entries(), ia, A.e[ia], t.a);
                value(v.m.data(), v.entries(), iv, v.e[iv], t.b);
                if(perturb_bits > 0 && ti == 0 && tb == 0)
                  {
                    mpf_div_2exp(t.d, t.a, (unsigned long)perturb_bits);
                    mpf_add(t.a, t.a, t.d);
                  }
                mpf_mul(t.c, t.a, t.b);
                mpf_add(sum, sum, t.c);
                mpf_abs(t.c, t.c);
                if(mpf_cmp(t.c, big) > 0)
                  mpf_set(big, t.c);
              }
          const int64_t out = b * g.ob + o * g.oo;
          value(O.m.data(), oent, out, O.e[out], t.a);
          mpf_sub(t.b, t.a, sum);
          const double err = log2_abs(t.b) - log2_abs(big);
          if(mpf_sgn(big) != 0 && err > worst)
            worst = err;
        }
    mpf_clears(sum, big, nullptr);
    report(name, worst, -double(NB) + 20);
  }
}

int main(int argc, char **argv)
{
  mpf_set_default_prec(16384);
  double perturb = 0;
  for(int i = 1; i + 1 < argc; i += 2)
    if(std::string(argv[i]) == "--perturb-bits")
      perturb = atof(argv[i + 1]);
    else if(std::string(argv[i]) == "--L")
      {
        L = atoi(argv[i + 1]);
        NB = 32 * L - 2;
      }

  Mpf_Vec t;
  // elementwise product
  {
    Host_Mat A = random_mat(2, 50, 40, 200, 0.1), B = random_mat(2, 50, 40, 200, 0.1);
    Dev dA(A), dB(B), dO(2, 50, 40);
    dev_mul(dA.b, dB.b, dO.b, 0);
    cudaStreamSynchronize(0);
    Host_Mat O = dO.get();
    double worst = -1e30;
    int bad_norm = 0;
    for(int64_t i = 0; i < A.entries(); ++i)
      {
        value(A.m.data(), A.entries(), i, A.e[i], t.a);
        value(B.m.data(), A.entries(), i, B.e[i], t.b);
        mpf_mul(t.c, t.a, t.b);
        value(O.m.data(), A.entries(), i, O.e[i], t.d);
        mpf_sub(t.a, t.d, t.c);
        if(mpf_sgn(t.c) != 0)
          worst = std::max(worst, log2_abs(t.a) - log2_abs(t.c));
        else if(mpf_sgn(t.d) != 0)
          ++bad_norm;
        if(O.e[i] != EZ)
          {
            uint32_t top = uint32_t(O.m[(L - 1) * A.entries() + i]);
            if(int32_t(top) < 0)
              top = ~top;
            if((top >> 29) != 1 && (top >> 29) != 0)
              ++bad_norm;
          }
      }
    report("mpf_mul, relative", worst, -double(NB) + 1);
    if(bad_norm)
      {
        printf("mpf_mul: %d not normalized\n", bad_norm);
        ++failures;
      }
  }
  // scalar product
  {
    Host_Mat A = random_mat(1, 30, 30, 50, 0.05), B = random_mat(1, 1, 1, 10, 0);
    Dev dA(A), dB(B), dO(1, 30, 30);
    dev_mul(dA.b, dB.b, dO.b, 0);
    cudaStreamSynchronize(0);
    Host_Mat O = dO.get();
    double worst = -1e30;
    value(B.m.data(), 1, 0, B.e[0], t.b);
    for(int64_t i = 0; i < A.entries(); ++i)
      {
        value(A.m.data(), A.entries(), i, A.e[i], t.a);
        mpf_mul(t.c, t.a, t.b);
        value(O.m.data(), A.entries(), i, O.e[i], t.d);
        mpf_sub(t.a, t.d, t.c);
        if(mpf_sgn(t.c) != 0)
          worst = std::max(worst, log2_abs(t.a) - log2_abs(t.c));
      }
    report("mpf_mul by a scalar, relative", worst, -double(NB) + 1);
  }
  // total, max_abs
  {
    Host_Mat A = random_mat(3, 40, 7, 400, 0.1);
    Dev dA(A);
    Host_Float s = mpf_total(dA.b, 0), mx = mpf_max_abs(dA.b, 0);
    mpf_t sum, big;
    mpf_inits(sum, big, nullptr);
    mpf_set_ui(sum, 0);
    mpf_set_ui(big, 0);
    for(int64_t i = 0; i < A.entries(); ++i)
      {
        value(A.m.data(), A.entries(), i, A.e[i], t.a);
        mpf_add(sum, sum, t.a);
        mpf_abs(t.a, t.a);
        if(mpf_cmp(t.a, big) > 0)
          mpf_set(big, t.a);
      }
    value(reinterpret_cast<int32_t *>(s.m.data()), 1, 0, s.e, t.b, int(s.m.size()));
    mpf_sub(t.c, t.b, sum);
    report("mpf_total, relative to the largest", log2_abs(t.c) - log2_abs(big),
           -double(NB) - 24);
    value(reinterpret_cast<int32_t *>(mx.m.data()), 1, 0, mx.e, t.b);
    mpf_sub(t.c, t.b, big);
    if(mpf_sgn(t.c) != 0)
      {
        printf("mpf_max_abs: differs from the exact maximum\n");
        ++failures;
      }
    else
      printf("mpf_max_abs                        exact  ok\n");
    mpf_clears(sum, big, nullptr);
  }
  // total: a negative term that rounds to zero is no term
  {
    Host_Mat A(1, 2, 1);
    set_random(A, 0, 0, 0.0);
    set_random(A, 1, 0, 0.0);
    // entry 0 positive, entry 1 negative
    const auto with_sign = [&](int64_t i, bool neg) {
      const bool is_neg = A.m[(L - 1) * A.entries() + i] < 0;
      if(is_neg != neg)
        {
          unsigned long long c = 1;
          for(int k = 0; k < L; ++k)
            {
              const unsigned long long v
                = (unsigned long long)(~uint32_t(A.m[k * A.entries() + i])) + c;
              A.m[k * A.entries() + i] = int32_t(uint32_t(v));
              c = v >> 32;
            }
        }
    };
    with_sign(0, false);
    with_sign(1, true);
    // t0 = e0 - 32, and a term is shifted by t0 - e <= 32 L: its top bit at NB - 1 is gone
    A.e[1] = A.e[0] - 32 - 32 * L;
    Dev dA(A);
    Host_Float s = mpf_total(dA.b, 0);
    value(A.m.data(), A.entries(), 0, A.e[0], t.a);
    value(reinterpret_cast<int32_t *>(s.m.data()), 1, 0, s.e, t.b, int(s.m.size()));
    mpf_sub(t.c, t.b, t.a);
    mpf_abs(t.d, t.a);
    report("mpf_total, a negative term rounding to zero", log2_abs(t.c) - log2_abs(t.d),
           -double(NB) - 24);
  }
  // sparse products against exact sums, and the slices of P
  {
    const int64_t n = 37, nv = 300;
    std::vector<int64_t> ptr(n + 1, 0);
    std::vector<uint32_t> col, val, out(n);
    for(int64_t o = 0; o < n; ++o)
      {
        const int terms = o == 5 ? 0 : 1 + int(rng() % 40);
        for(int k = 0; k < terms; ++k)
          {
            col.push_back(uint32_t(rng() % nv) | ((rng() & 1) ? 0x80000000u : 0u));
            val.push_back(0);
          }
        ptr[o + 1] = int64_t(col.size());
        out[o] = uint32_t(2 * o + 1);
      }
    const int64_t nnz = int64_t(col.size());
    for(int64_t k = 0; k < nnz; ++k)
      val[k] = uint32_t((nnz - 1) - k);
    Host_Mat A = random_mat(1, nnz, 1, 300, 0.05), v = random_mat(1, nv, 1, 300, 0.05);
    Dev dA(A), dv(v), dO(1, 2 * n + 1, 1);
    mpf_zero(dO.b, 0);
    auto dev = [](const void *p, size_t bytes) {
      void *d = nullptr;
      if(cudaMalloc(&d, bytes) != cudaSuccess || cudaMemcpy(d, p, bytes, cudaMemcpyHostToDevice) != cudaSuccess)
        throw std::runtime_error("device copy");
      return d;
    };
    Csr_Gemv g;
    g.n = n;
    g.ptr = static_cast<const int64_t *>(dev(ptr.data(), ptr.size() * 8));
    g.col = static_cast<const uint32_t *>(dev(col.data(), col.size() * 4));
    g.val = static_cast<const uint32_t *>(dev(val.data(), val.size() * 4));
    g.out = static_cast<const uint32_t *>(dev(out.data(), out.size() * 4));
    mpf_csr_gemv(dA.b, dv.b, dO.b, g, 0);
    cudaStreamSynchronize(0);
    Host_Mat O = dO.get();
    double worst = -1e30, worst_control = -1e30;
    mpf_t sum, big;
    mpf_inits(sum, big, nullptr);
    for(int64_t o = 0; o < n; ++o)
      {
        mpf_set_ui(sum, 0);
        mpf_set_ui(big, 0);
        for(int64_t k = ptr[o]; k < ptr[o + 1]; ++k)
          {
            value(A.m.data(), A.entries(), val[k], A.e[val[k]], t.a);
            const int64_t iv = col[k] & 0x7fffffffu;
            value(v.m.data(), v.entries(), iv, v.e[iv], t.b);
            mpf_mul(t.c, t.a, t.b);
            if(col[k] & 0x80000000u)
              mpf_mul_2exp(t.c, t.c, 1);
            mpf_add(sum, sum, t.c);
            mpf_abs(t.c, t.c);
            if(mpf_cmp(t.c, big) > 0)
              mpf_set(big, t.c);
          }
        const int64_t at = out[o];
        value(O.m.data(), O.entries(), at, O.e[at], t.a);
        mpf_sub(t.b, t.a, sum);
        if(mpf_sgn(big) != 0)
          worst = std::max(worst, log2_abs(t.b) - log2_abs(big));
        // the even outputs are not written
        value(O.m.data(), O.entries(), 2 * o, O.e[2 * o], t.a);
        worst_control = std::max(worst_control, log2_abs(t.a));
      }
    mpf_clears(sum, big, nullptr);
    report("mpf_csr_gemv, relative to the largest term", worst, -double(NB) + 20);
    if(worst_control > -1e29)
      {
        printf("mpf_csr_gemv: wrote outside its output list\n");
        ++failures;
      }
    // slices: entries of two 4 x 4 blocks of three slices, with their mirrors, and P's rows
    const int64_t bn = 4;
    const int dp = 3, nb = 2;
    std::vector<uint32_t> pos, slice;
    std::vector<int32_t> vm;
    Host_Mat F = random_mat(1, 5, 1, 10, 0.0);
    const int64_t ents[5][4] = {{0, 0, 2, 1}, {0, 2, 3, 3}, {1, 1, 1, 0}, {1, 0, 0, 0}, {0, 1, 3, 0}}; // z, s, i, j
    for(auto &e : ents)
      {
        pos.push_back(uint32_t(e[0] * bn * bn + e[3] * bn + e[2]) | 0x80000000u);
        slice.push_back(uint32_t(7 + e[1])); // p0 = 7
      }
    Dev dF(F), dW(nb * dp, bn, bn);
    mpf_zero(dW.b, 0);
    const auto *dpos = static_cast<const uint32_t *>(dev(pos.data(), pos.size() * 4));
    const auto *dslice = static_cast<const uint32_t *>(dev(slice.data(), slice.size() * 4));
    mpf_fill_slices(dF.b, dpos, dslice, 0, 5, bn, dp, 7, dW.b, 0);
    cudaStreamSynchronize(0);
    Host_Mat W = dW.get();
    int bad = 0;
    for(int k = 0; k < 5; ++k)
      {
        const int64_t base = (ents[k][0] * dp + ents[k][1]) * bn * bn;
        for(const int64_t q : {base + ents[k][3] * bn + ents[k][2], base + ents[k][2] * bn + ents[k][3]})
          {
            if(W.e[q] != F.e[k])
              ++bad;
            for(int w = 0; w < L; ++w)
              if(W.m[w * W.entries() + q] != F.m[w * F.entries() + k])
                ++bad;
          }
      }
    int nonzero = 0;
    for(int64_t q = 0; q < W.entries(); ++q)
      nonzero += W.e[q] != EZ;
    if(nonzero != 8) // five entries, two on the diagonal
      ++bad;
    // the valid rows of block 1 (dimension 2) and of block 0 (dimension 3)
    std::vector<uint32_t> tab;
    for(int z = 0; z < nb; ++z)
      for(int j = 0; j < 3 - z; ++j)
        for(int i = 0; i < 3 - z; ++i)
          tab.push_back(uint32_t(z * bn * bn + j * bn + i));
    const int64_t R = int64_t(tab.size());
    const auto *dtab = static_cast<const uint32_t *>(dev(tab.data(), tab.size() * 4));
    Dev dP(1, R, 10);
    mpf_zero(dP.b, 0);
    mpf_gather_slices(dW.b, dtab, R, bn, dp, 7, dP.b, 0);
    cudaStreamSynchronize(0);
    Host_Mat Pm = dP.get();
    for(int s_ = 0; s_ < dp; ++s_)
      for(int64_t r = 0; r < R; ++r)
        {
          const int64_t z = tab[r] / (bn * bn), rem = tab[r] % (bn * bn);
          const int64_t from = (z * dp + s_) * bn * bn + rem, to = (7 + s_) * R + r;
          if(Pm.e[to] != W.e[from])
            ++bad;
          for(int w = 0; w < L; ++w)
            if(Pm.m[w * Pm.entries() + to] != W.m[w * W.entries() + from])
              ++bad;
        }
    printf("mpf_fill_slices, mpf_gather_slices: %d differences  %s\n", bad, bad ? "FAILED" : "ok");
    if(bad)
      ++failures;
  }
  // to double
  {
    Host_Mat A = random_mat(3, 20, 20, 300, 0.1);
    Dev dA(A);
    double *d;
    cudaMalloc(&d, A.entries() * 8);
    std::vector<int64_t> scale(3);
    mpf_to_double(dA.b, d, scale.data(), 0);
    cudaStreamSynchronize(0);
    std::vector<double> h(A.entries());
    cudaMemcpy(h.data(), d, h.size() * 8, cudaMemcpyDeviceToHost);
    cudaFree(d);
    double worst = -1e30;
    for(int64_t i = 0; i < A.entries(); ++i)
      {
        const int b = int(i / 400);
        value(A.m.data(), A.entries(), i, A.e[i], t.a);
        mpf_set_d(t.b, h[i]);
        if(scale[b] >= 0)
          mpf_mul_2exp(t.b, t.b, scale[b]);
        else
          mpf_div_2exp(t.b, t.b, -scale[b]);
        mpf_sub(t.c, t.b, t.a);
        // the double clamps entries below 2^-1100 of the block's largest
        if(mpf_sgn(t.a) != 0 && log2_abs(t.a) - double(scale[b]) > -1000)
          worst = std::max(worst, log2_abs(t.c) - log2_abs(t.a));
      }
    report("mpf_to_double, relative", worst, -52);
  }
  // matrix-vector products
  {
    const int nb = 3;
    const int64_t P = 37, N = 29;
    Host_Mat B = random_mat(nb, P, N, 30, 0.1), y = random_mat(1, N, 1, 30, 0.1);
    Gemv_Spec g;
    g.batch = nb;
    g.no = P;
    g.ab = P * N;
    g.ao = 1;
    g.ati = P;
    g.nti = N;
    g.vti = 1;
    g.ob = P;
    g.oo = 1;
    check_gemv("gemv B y (batched)", B, y, g, perturb);

    Host_Mat x = random_mat(nb, P, 1, 30, 0.1);
    Gemv_Spec h;
    h.no = N;
    h.ao = P;
    h.atb = P * N;
    h.ati = 1;
    h.ntb = nb;
    h.nti = P;
    h.vtb = P;
    h.vti = 1;
    h.ob = N;
    h.oo = 1;
    check_gemv("gemv B^T x (over blocks)", B, x, h, perturb);

    Host_Mat M = random_mat(nb, P, P, 300, 0.2), v = random_mat(nb, P, 1, 300, 0.1);
    Gemv_Spec s;
    s.batch = nb;
    s.no = P;
    s.ab = P * P;
    s.ao = 1;
    s.ati = P;
    s.nti = P;
    s.vb = P;
    s.vti = 1;
    s.ob = P;
    s.oo = 1;
    check_gemv("gemv M v, wide exponents", M, v, s, perturb);
    Gemv_Spec st = s;
    st.ao = P;
    st.ati = 1;
    check_gemv("gemv M^T v, wide exponents", M, v, st, perturb);
  }
  // B as its nonzeros: unpacking gives back the same bytes, including the entries
  // that only look like zeros; a wrong position must fail this
  {
    Host_Mat A = random_mat(3, 40, 50, 200, 0.8);
    const int64_t tot = A.entries();
    A.e[5] = MPF_EZ;       // exponent of a zero, mantissa not
    A.m[size_t(3) * tot + 5] = 1;
    A.e[6] = 17;           // mantissa of a zero, exponent not
    for(int t = 0; t < L; ++t)
      A.m[size_t(t) * tot + 6] = 0;
    std::vector<uint32_t> pos;
    std::vector<int32_t> vm;
    std::vector<int64_t> ve;
    mpf_pack_sparse(A.m.data(), A.e.data(), tot, L, pos, vm, ve);
    const int64_t nnz = int64_t(pos.size());
    Host_Mat Hv(1, nnz, 1);
    Hv.m = vm;
    Hv.e = ve;
    auto roundtrip = [&](const std::vector<uint32_t> &p) {
      Dev dv(Hv), dO(A.batch, A.rows, A.cols);
      uint32_t *dp = nullptr;
      if(cudaMalloc(&dp, p.size() * 4) != cudaSuccess)
        throw std::runtime_error("cudaMalloc");
      cudaMemcpy(dp, p.data(), p.size() * 4, cudaMemcpyHostToDevice);
      mpf_unpack_sparse(dv.b, dp, dO.b, 0);
      cudaStreamSynchronize(0);
      Host_Mat O = dO.get();
      cudaFree(dp);
      return O.m == A.m && O.e == A.e;
    };
    std::vector<uint32_t> wrong = pos;
    wrong[nnz / 2] = wrong[nnz / 2] + 1 == wrong[nnz / 2 + 1]
                       ? wrong[nnz / 2] - 1
                       : wrong[nnz / 2] + 1;
    const bool same = roundtrip(pos), control = roundtrip(wrong);
    const bool ok = same && !control && nnz > tot / 10 && nnz < tot / 2;
    printf("sparse B round trip, %ld of %ld entries kept: %s, wrong position %s  %s\n",
           long(nnz), long(tot), same ? "same bytes" : "DIFFERS",
           control ? "NOT CAUGHT" : "caught", ok ? "ok" : "FAILED");
    if(!ok)
      ++failures;
  }
  printf("%s\n", failures ? "FAILED" : "all ok");
  return failures ? 1 : 0;
}
