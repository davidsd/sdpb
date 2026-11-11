// GPU backend for SDPB: check and benchmark of the batched float product
// Mpf_Gemm, C_b = A_b^T B_b (floats in, floats out).  Random batched floats
// (L limbs, column scales spread over 2^-32..2^31, 1/16 zeros); sampled (or
// all) entries compared with the exact sums computed by GMP from the same
// mantissas and exponents; error relative to k max|A_i| max|B_j|.
//
//   mpf_gemm_check --m 945 --n 3406 --k 945 --batch 36 --L 28 --digits 112
//                  [--panel 512] [--samples 20] [--reps 2] [--timing]

#include "../rns/Mpf_Gram.hxx"

#include <gmp.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace sdpb_gpu;

namespace
{
  __host__ __device__ inline uint64_t splitmix(uint64_t x)
  {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
  }

  // entry (i, j) of chunk c: limbs into w[0..L), returns the exponent
  __host__ __device__ inline int64_t gen_entry(uint64_t seed, int64_t c,
                                               int64_t i, int64_t j, int L,
                                               uint32_t *w, int triangle = 0)
  {
    uint64_t h = splitmix(seed ^ splitmix(uint64_t(c) * 0x100000001B3ULL
                                          ^ splitmix(uint64_t(j) << 22
                                                     ^ uint64_t(i))));
    if((h & 15) == 0 || (triangle == 2 && i < j)
       || ((triangle == 1 || triangle == 4) && i > j))
      {
        for(int t = 0; t < L; ++t)
          w[t] = 0;
        return MPF_EZ;
      }
    for(int t = 0; t < L - 1; ++t)
      w[t] = uint32_t(splitmix(h + 1 + t));
    w[L - 1] = (uint32_t(splitmix(h + 1 + L)) & 0x1FFFFFFFu) | 0x20000000u;
    if((h >> 8) & 1)
      {
        unsigned long long cc = 1;
        for(int t = 0; t < L; ++t)
          {
            unsigned long long v = (unsigned long long)(~w[t]) + cc;
            w[t] = uint32_t(v);
            cc = v >> 32;
          }
      }
    int64_t col_scale = int64_t((j * 37) % 64) - 32;
    return -(32LL * L - 2) + int64_t((h >> 16) % 9) - 4 + col_scale;
  }

  // batched: block b of (rows x cols), entry (i, j) -> gen_entry(seed, b, i, j)
  __global__ void fill(int32_t *m, int64_t *e, int64_t rows, int64_t N, int L,
                       uint64_t seed, int64_t batch, int triangle)
  {
    int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    if(idx >= batch * rows * N)
      return;
    int64_t b = idx / (rows * N), rem = idx % (rows * N);
    uint32_t w[64];
    e[idx] = gen_entry(seed, b, rem % rows, rem / rows, L, w, triangle);
    for(int t = 0; t < L; ++t)
      m[(int64_t)t * batch * rows * N + idx] = int32_t(w[t]);
  }

  void tc_to_mpz(mpz_t x, const uint32_t *w, int n)
  {
    mpz_import(x, n, -1, 4, 0, 0, w);
    if(w[n - 1] >> 31)
      {
        mpz_t t;
        mpz_init(t);
        mpz_ui_pow_ui(t, 2, 32 * n);
        mpz_sub(x, x, t);
        mpz_clear(t);
      }
  }

  // log2 |x * 2^e|
  double log2_abs(const mpz_t x, int64_t e)
  {
    if(mpz_sgn(x) == 0)
      return -1e300;
    long ex = 0;
    double d = mpz_get_d_2exp(&ex, x);
    return std::log2(std::fabs(d)) + double(ex) + double(e);
  }

  const char *arg(int argc, char **argv, const char *name, const char *dflt)
  {
    for(int i = 1; i + 1 < argc; ++i)
      if(!strcmp(argv[i], name))
        return argv[i + 1];
    return dflt;
  }
  bool flag(int argc, char **argv, const char *name)
  {
    for(int i = 1; i < argc; ++i)
      if(!strcmp(argv[i], name))
        return true;
    return false;
  }
}

int main(int argc, char **argv)
{
  const int m = atoi(arg(argc, argv, "--m", "40"));
  const int n = atoi(arg(argc, argv, "--n", "70"));
  const int k = atoi(arg(argc, argv, "--k", "50"));
  const int batch = atoi(arg(argc, argv, "--batch", "3"));
  const int L = atoi(arg(argc, argv, "--L", "28"));
  const int D = atoi(arg(argc, argv, "--digits", "112"));
  const int panel = atoi(arg(argc, argv, "--panel", "512"));
  const int samples = atoi(arg(argc, argv, "--samples", "0"));
  const int reps = atoi(arg(argc, argv, "--reps", "1"));
  const int triangles = atoi(arg(argc, argv, "--triangles", "0"));
  const int triA = triangles == 1 ? 1 : 0, triB = triangles >= 2 ? triangles : 0;
  Gram_Options opt;
  opt.prime_group = atoi(arg(argc, argv, "--group", "32"));
  opt.timing = flag(argc, argv, "--timing");
  opt.crt_base128 = !flag(argc, argv, "--crt256");
  opt.fast_rounding = !flag(argc, argv, "--slow-rounding");
  opt.reduce_to_crt = !flag(argc, argv, "--no-reduce-crt");
  // --compare-rounding: the product again with the other rounding kernel,
  // every limb and exponent compared
  const bool cmp_round = flag(argc, argv, "--compare-rounding");
  opt.crt_batch = atoll(arg(argc, argv, "--crt-batch", "524288"));
  opt.fused_residues = true;
  // --crt-guard g: the CRT forms only the top of each product. The gate below
  // must hold for g > 0; a negative g cuts into the float's own bits, and the
  // gate must fail (the control)
  opt.crt_guard_bits = atoi(arg(argc, argv, "--crt-guard", "0"));
  const uint64_t sA = 101, sB = 202;
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("device %s; Mpf_Gemm batch=%d of (%d x %d)^T (%d x %d), L=%d, digits %d, "
         "panel %d\n", prop.name, batch, k, m, k, n, L, D, panel);
  Mpf_Matrix A, B;
  A.L = B.L = L;
  A.rows = B.rows = k;
  A.cols = m;
  B.cols = n;
  const int64_t eA = int64_t(batch) * k * m, eB = int64_t(batch) * k * n,
                eC = int64_t(batch) * m * n;
  int32_t *cm;
  int64_t *ce;
  if(cudaMalloc(&A.m, eA * L * 4) || cudaMalloc(&A.e, eA * 8)
     || cudaMalloc(&B.m, eB * L * 4) || cudaMalloc(&B.e, eB * 8)
     || cudaMalloc(&cm, eC * L * 4) || cudaMalloc(&ce, eC * 8))
    {
      fprintf(stderr, "cannot allocate\n");
      return 2;
    }
  fill<<<(eA + 255) / 256, 256>>>(A.m, A.e, k, m, L, sA, batch, triA);
  fill<<<(eB + 255) / 256, 256>>>(B.m, B.e, k, n, L, sB, batch, triB);
  cudaDeviceSynchronize();
  Mpf_Gemm G(m, n, k, batch, L, D, opt, panel);
  printf("primes %zu, device buffers %.2f GB\n", G.gemm().basis().primes.size(),
         G.gemm().device_bytes() / 1e9);
  double best = 1e30;
  for(int rep = 0; rep < reps; ++rep)
    {
      cudaDeviceSynchronize();
      auto t0 = std::chrono::steady_clock::now();
      G.multiply(A, B, cm, ce, triangles);
      cudaDeviceSynchronize();
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      printf("rep %d: %.4f s (floats in, floats out)\n", rep, s);
      best = std::min(best, s);
    }
  double mpmac = double(batch) * m * n * k;
  printf("best %.4f s; %.3g MP-MAC -> %.3g MP-MAC/s\n", best, mpmac, mpmac / best);
  for(auto &[kk, v] : G.timings())
    printf("  phase %-22s %9.2f ms (all reps)\n", kk.c_str(), v);
  for(auto &[kk, v] : G.gemm().timings())
    printf("  phase %-22s %9.2f ms (all reps)\n", kk.c_str(), v);
  // check sampled entries against exact GMP sums
  std::vector<int32_t> hm(size_t(L) * eC);
  std::vector<int64_t> he(eC);
  cudaMemcpy(hm.data(), cm, hm.size() * 4, cudaMemcpyDeviceToHost);
  cudaMemcpy(he.data(), ce, he.size() * 8, cudaMemcpyDeviceToHost);
  bool round_ok = true;
  if(cmp_round || triangles)
    {
      Gram_Options o2 = opt;
      if(cmp_round)
        o2.fast_rounding = !opt.fast_rounding;
      o2.timing = false;
      {
        Mpf_Gemm G2(m, n, k, batch, L, D, o2, panel);
        G2.multiply(A, B, cm, ce);
        cudaDeviceSynchronize();
      }
      std::vector<int32_t> hm2(hm.size());
      std::vector<int64_t> he2(he.size());
      cudaMemcpy(hm2.data(), cm, hm2.size() * 4, cudaMemcpyDeviceToHost);
      cudaMemcpy(he2.data(), ce, he2.size() * 8, cudaMemcpyDeviceToHost);
      int64_t diff = 0;
      for(int64_t e = 0; e < eC; ++e)
        {
          bool same = he[e] == he2[e];
          for(int t = 0; t < L && same; ++t)
            same = hm[size_t(t) * eC + e] == hm2[size_t(t) * eC + e];
          diff += !same;
        }
      round_ok = diff == 0;
      printf("general/alternate product compared: %lld of %lld entries differ\n",
             (long long)diff, (long long)eC);
    }
  std::mt19937_64 rng(9);
  std::vector<int64_t> idx;
  if(samples == 0)
    for(int64_t e = 0; e < eC; ++e)
      idx.push_back(e);
  else
    {
      idx = {0, eC - 1, int64_t(m) - 1, int64_t(m) * n};
      while((int)idx.size() < samples)
        idx.push_back(int64_t(rng() % eC));
    }
  std::vector<uint32_t> w(64), w2(64);
  mpz_t S, t, O, a, bz;
  mpz_inits(S, t, O, a, bz, nullptr);
  double worst = -1e300;
  long bad = 0;
  const double nbits = std::min(8.0 * D - 2, 32.0 * L - 2);
  for(int64_t eo : idx)
    {
      if(eo < 0 || eo >= eC)
        continue;
      int64_t b = eo / (int64_t(m) * n), rem = eo % (int64_t(m) * n);
      int64_t j = rem / m, i = rem % m;
      // column maxima of A(:, i) and B(:, j), and the exact sum
      double cA = -1e300, cB = -1e300;
      std::vector<int64_t> ea(k), eb(k);
      std::vector<__mpz_struct> za(k), zb(k);
      int64_t Emin = INT64_MAX;
      for(int r = 0; r < k; ++r)
        {
          mpz_init(&za[r]);
          mpz_init(&zb[r]);
          ea[r] = gen_entry(sA, b, r, i, L, w.data(), triA);
          if(ea[r] != MPF_EZ) { tc_to_mpz(&za[r], w.data(), L); cA = std::max(cA, log2_abs(&za[r], ea[r])); }
          eb[r] = gen_entry(sB, b, r, j, L, w2.data(), triB);
          if(eb[r] != MPF_EZ) { tc_to_mpz(&zb[r], w2.data(), L); cB = std::max(cB, log2_abs(&zb[r], eb[r])); }
          if(ea[r] != MPF_EZ && eb[r] != MPF_EZ)
            Emin = std::min(Emin, ea[r] + eb[r]);
        }
      mpz_set_ui(S, 0);
      for(int r = 0; r < k && Emin != INT64_MAX; ++r)
        if(ea[r] != MPF_EZ && eb[r] != MPF_EZ)
          {
            mpz_mul(t, &za[r], &zb[r]);
            mpz_mul_2exp(t, t, ea[r] + eb[r] - Emin);
            mpz_add(S, S, t);
          }
      for(int r = 0; r < k; ++r) { mpz_clear(&za[r]); mpz_clear(&zb[r]); }
      for(int q = 0; q < L; ++q)
        w[q] = uint32_t(hm[size_t(q) * eC + eo]);
      int64_t oe = he[eo];
      if(oe == MPF_EZ) mpz_set_ui(O, 0); else tc_to_mpz(O, w.data(), L);
      if(Emin == INT64_MAX) Emin = oe == MPF_EZ ? 0 : oe;
      int64_t ez = oe == MPF_EZ ? Emin : std::min(oe, Emin);
      mpz_mul_2exp(t, O, oe == MPF_EZ ? 0 : oe - ez);
      mpz_mul_2exp(S, S, Emin - ez);
      mpz_sub(t, t, S);
      double rn = mpz_sgn(t) == 0 ? -1e300
                    : log2_abs(t, ez) - (cA + cB + std::log2(double(k)));
      worst = std::max(worst, rn);
      if(rn > -(nbits - 4)) ++bad;
    }
  if(!round_ok)
    ++bad;
  printf("CHECK %s: %zu entries; worst error / (k max|A_i| max|B_j|) = 2^%.1f (gate 2^-%.0f)%s\n",
         bad ? "FAILED" : "PASSED", idx.size(), worst, nbits - 4,
         round_ok ? "" : "; rounding kernels differ");
  return bad ? 1 : 0;
}
