// GPU backend for SDPB: exactness check and benchmark of Rns_Gemm (batched
// exact C_b = A_b^T B_b).  A, B from a counter-based hash (balanced digits),
// on the device and again on the host; sampled entries compared with GMP.
//
//   gemm_check --m 945 --n 3406 --k 944 --batch 36 --digits 112 [--panel 512]
//              [--samples 200] [--reps 2] [--timing]
// (the nmax14 dim-5 Schur blocks: P_j = 945, k padded to a multiple of 16)

#include "../rns/Rns_Gemm.hxx"

#include <gmp.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
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
  __host__ __device__ inline int8_t digit(uint64_t seed, int64_t entry, int t)
  {
    return (int8_t)(splitmix(splitmix(seed ^ splitmix(uint64_t(entry))) + t)
                    & 255);
  }
  __global__ void fill(int8_t *d, int64_t entries, int D, uint64_t seed)
  {
    int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    if(e >= entries)
      return;
    for(int t = 0; t < D; ++t)
      d[e * D + t] = digit(seed, e, t);
  }
  void entry_mpz(mpz_t a, uint64_t seed, int64_t e, int D)
  {
    mpz_set_si(a, digit(seed, e, D - 1));
    for(int t = D - 2; t >= 0; --t)
      {
        mpz_mul_2exp(a, a, 8);
        long d = digit(seed, e, t);
        if(d >= 0)
          mpz_add_ui(a, a, d);
        else
          mpz_sub_ui(a, a, -d);
      }
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
  const int m = atoi(arg(argc, argv, "--m", "64"));
  const int n = atoi(arg(argc, argv, "--n", "100"));
  const int k = atoi(arg(argc, argv, "--k", "64"));
  const int batch = atoi(arg(argc, argv, "--batch", "3"));
  const int DA = atoi(arg(argc, argv, "--digits", "112"));
  const int DB = atoi(arg(argc, argv, "--digitsB", arg(argc, argv, "--digits", "112")));
  const int panel = atoi(arg(argc, argv, "--panel", "512"));
  const int samples = atoi(arg(argc, argv, "--samples", "0"));
  const int reps = atoi(arg(argc, argv, "--reps", "1"));
  Gram_Options opt;
  opt.prime_group = atoi(arg(argc, argv, "--group", "16"));
  opt.timing = flag(argc, argv, "--timing");
  opt.crt_batch = atoll(arg(argc, argv, "--crt-batch", "524288"));
  opt.fused_residues = flag(argc, argv, "--fused");
  opt.fused_crt = flag(argc, argv, "--fused-crt");
  opt.crt_base128 = !flag(argc, argv, "--crt256");
  opt.reduce_to_crt = !flag(argc, argv, "--no-reduce-crt");
  const uint64_t sA = 11, sB = 23;
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("device %s; Rns_Gemm batch=%d of (%d x %d)^T (%d x %d), digits %d/%d, "
         "panel %d\n",
         prop.name, batch, k, m, k, n, DA, DB, panel);
  const int64_t entA = int64_t(batch) * k * m, entB = int64_t(batch) * k * n;
  int8_t *dA, *dB;
  if(cudaMalloc(&dA, entA * DA) || cudaMalloc(&dB, entB * DB))
    {
      fprintf(stderr, "cannot allocate inputs\n");
      return 2;
    }
  fill<<<(entA + 255) / 256, 256>>>(dA, entA, DA, sA);
  fill<<<(entB + 255) / 256, 256>>>(dB, entB, DB, sB);
  Rns_Gemm G(m, n, k, batch, DA, DB, opt, panel);
  const int W = G.limbs_out();
  int32_t *dC;
  if(cudaMalloc(&dC, size_t(W) * G.plane() * 4))
    {
      fprintf(stderr, "cannot allocate output\n");
      return 2;
    }
  printf("primes %zu, output limbs %d, device buffers %.2f GB, inputs %.2f GB, "
         "output %.2f GB\n",
         G.basis().primes.size(), W, G.device_bytes() / 1e9,
         (entA * DA + entB * DB) / 1e9, W * 4.0 * G.plane() / 1e9);
  double best = 1e30;
  for(int rep = 0; rep < reps; ++rep)
    {
      cudaDeviceSynchronize();
      double ops0 = G.int8_ops();
      auto t0 = std::chrono::steady_clock::now();
      G.multiply(dA, dB, dC);
      cudaDeviceSynchronize();
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                               - t0)
                   .count();
      printf("rep %d: %.4f s, %.0f TOPS effective (all INT8 products)\n", rep,
             s, (G.int8_ops() - ops0) / s / 1e12);
      best = std::min(best, s);
    }
  double mpmac = double(batch) * m * n * k;
  printf("best %.4f s; %.3g MP-MAC -> %.3g MP-MAC/s\n", best, mpmac,
         mpmac / best);
  for(auto &[kk, v] : G.timings())
    printf("  phase %-12s %9.2f ms (all reps)\n", kk.c_str(), v);

  // check
  std::vector<int32_t> h(size_t(W) * G.plane());
  cudaMemcpy(h.data(), dC, h.size() * 4, cudaMemcpyDeviceToHost);
  std::mt19937_64 rng(7);
  const int64_t total = G.plane();
  std::vector<int64_t> idx;
  if(samples == 0)
    for(int64_t e = 0; e < total; ++e)
      idx.push_back(e);
  else
    {
      idx = {0, total - 1, int64_t(m) - 1, int64_t(m) * n - 1,
             int64_t(m) * n};
      while((int)idx.size() < samples)
        idx.push_back(int64_t(rng() % total));
    }
  mpz_t a, b, ref, got, two;
  mpz_inits(a, b, ref, got, two, nullptr);
  mpz_ui_pow_ui(two, 2, 32 * W);
  long bad = 0;
  std::vector<uint32_t> w(W);
  for(int64_t e : idx)
    {
      if(e < 0 || e >= total)
        continue;
      int64_t bb = e / (int64_t(m) * n), rem = e % (int64_t(m) * n);
      int64_t j = rem / m, i = rem % m;
      mpz_set_ui(ref, 0);
      for(int t = 0; t < k; ++t)
        {
          entry_mpz(a, sA, (bb * m + i) * k + t, DA);
          entry_mpz(b, sB, (bb * n + j) * k + t, DB);
          mpz_addmul(ref, a, b);
        }
      for(int q = 0; q < W; ++q)
        w[q] = uint32_t(h[size_t(q) * total + e]);
      mpz_import(got, W, -1, 4, 0, 0, w.data());
      if(w[W - 1] >> 31)
        mpz_sub(got, got, two);
      if(mpz_cmp(ref, got))
        {
          if(bad < 3)
            gmp_printf("MISMATCH block %ld (%ld,%ld): ref %Zd got %Zd\n",
                       (long)bb, (long)i, (long)j, ref, got);
          ++bad;
        }
    }
  printf("CHECK %s: %ld of %zu entries differ from GMP\n",
         bad ? "FAILED" : "PASSED", bad, idx.size());
  return bad ? 1 : 0;
}
