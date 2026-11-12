// GPU backend for SDPB: check and benchmark of the Q-formation path on device
// floats (Mpf_Gram): random P chunks (rows x N, L-limb floats with column scales
// spread over 2^-32..2^31 and 1/16 zeros) -> Q = sum P^T P (lower triangle,
// floats).  Sampled entries are compared with the exact sum computed by GMP
// (integer mantissas and exponents, independent of the residue path); the
// error is reported relative to the normwise bound rows * max|P_i| * max|P_j|.
//
//   mpf_gram_check --N 3406 --rows 17520 --chunks 4 --L 28 --digits 112
//                  --sample-cols 12 [--reps 2] [--timing]

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
                                               uint32_t *w)
  {
    uint64_t h = splitmix(seed ^ splitmix(uint64_t(c) * 0x100000001B3ULL
                                          ^ splitmix(uint64_t(j) << 22
                                                     ^ uint64_t(i))));
    // --bands B (in the seed's top byte): rows in runs of 64 share a band b =
    // (i / 64) % B, and their entries in the columns with j % B > b are exact
    // zeros, as the rows of SDPB's P = L^-1 B are zero in whole groups of
    // columns. Then no random zeros, so that the rows of a band share their
    // nonzero columns.
    const int B = int(seed >> 56);
    if(B ? (j % B) > ((i / 64) % B) : (h & 15) == 0)
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

  __global__ void fill(int32_t *m, int64_t *e, int64_t rows, int64_t N, int L,
                       uint64_t seed, int64_t c)
  {
    int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    if(idx >= rows * N)
      return;
    uint32_t w[64];
    e[idx] = gen_entry(seed, c, idx % rows, idx / rows, L, w); // column-major
    for(int t = 0; t < L; ++t)
      m[(int64_t)t * rows * N + idx] = int32_t(w[t]);
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
  const int64_t N = atoll(arg(argc, argv, "--N", "200"));
  const int64_t rows = atoll(arg(argc, argv, "--rows", "1000"));
  const int chunks = atoi(arg(argc, argv, "--chunks", "2"));
  const int L = atoi(arg(argc, argv, "--L", "28"));
  const int D = atoi(arg(argc, argv, "--digits", "112"));
  const int sample_cols = atoi(arg(argc, argv, "--sample-cols", "0"));
  const int reps = atoi(arg(argc, argv, "--reps", "1"));
  // --structured: Q by groups of rows with the same nonzero columns
  // (Mpf_Gram::add_structured over all chunks) instead of add(); --plane: its
  // bound on a residue plane's entries
  const bool structured = flag(argc, argv, "--structured");
  const int64_t plane = atoll(arg(argc, argv, "--plane", "0"));
  const int bands = atoi(arg(argc, argv, "--bands", "0"));
  // negative control of --structured: one group of rows is left out
  const int skip_group = atoi(arg(argc, argv, "--skip-group", "-1"));
  const uint64_t seed = (strtoull(arg(argc, argv, "--seed", "5"), nullptr, 10)
                         & ((uint64_t(1) << 56) - 1))
                        | (uint64_t(bands) << 56);
  // --reuse: every chunk is chunk 0 (one buffer)
  const bool reuse = flag(argc, argv, "--reuse");
  Gram_Options opt;
  opt.prime_group = atoi(arg(argc, argv, "--group", "16"));
  opt.timing = flag(argc, argv, "--timing");
  opt.crt_batch = atoll(arg(argc, argv, "--crt-batch", "524288"));
  opt.fused_residues = flag(argc, argv, "--fused");
  opt.fused_crt = flag(argc, argv, "--fused-crt");

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("device %s; Mpf_Gram N=%lld rows/chunk=%lld chunks=%d L=%d (%d-bit "
         "mantissas) digits=%d%s\n",
         prop.name, (long long)N, (long long)rows, chunks, L, 32 * L - 2, D,
         reuse ? " (one chunk reused)" : "");
  const int nbuf = reuse ? 1 : chunks;
  std::vector<Mpf_Matrix> Pc(nbuf);
  for(int c = 0; c < nbuf; ++c)
    {
      Pc[c].L = L;
      Pc[c].rows = rows;
      Pc[c].cols = N;
      if(cudaMalloc(&Pc[c].m, size_t(L) * rows * N * 4) != cudaSuccess
         || cudaMalloc(&Pc[c].e, size_t(rows) * N * 8) != cudaSuccess)
        {
          fprintf(stderr, "cannot allocate chunk %d\n", c);
          return 2;
        }
      fill<<<(rows * N + 255) / 256, 256>>>(Pc[c].m, Pc[c].e, rows, N, L,
                                            seed, c);
    }
  cudaDeviceSynchronize();
  Mpf_Gram MG(N, L, D, rows * chunks, opt);
  const int64_t total = MG.gram().packed_size();
  int32_t *d_qm = nullptr;
  int64_t *d_qe = nullptr;
  if(cudaMalloc(&d_qm, size_t(L) * total * 4) != cudaSuccess
     || cudaMalloc(&d_qe, size_t(total) * 8) != cudaSuccess)
    {
      fprintf(stderr, "cannot allocate Q\n");
      return 2;
    }
  printf("primes %zu, device buffers %.2f GB (+ chunks %.2f GB, Q %.2f GB)\n",
         MG.gram().basis().primes.size(), MG.gram().device_bytes() / 1e9,
         nbuf * (L * 4.0 + 8) * rows * N / 1e9, (L * 4.0 + 8) * total / 1e9);
  double best = 1e30;
  for(int rep = 0; rep < reps; ++rep)
    {
      cudaDeviceSynchronize();
      auto t0 = std::chrono::steady_clock::now();
      for(int c = 0; c < chunks; ++c)
        MG.scan(Pc[reuse ? 0 : c]);
      if(structured)
        {
          std::vector<Mpf_Matrix> all;
          for(int c = 0; c < chunks; ++c)
            all.push_back(Pc[reuse ? 0 : c]);
          MG.add_structured(all, plane, skip_group);
          printf("structured: %.4f of the dense product's multiply-adds\n",
                 MG.structured_fraction());
        }
      else
        for(int c = 0; c < chunks; ++c)
          MG.add(Pc[reuse ? 0 : c]);
      MG.finish(d_qm, d_qe);
      cudaDeviceSynchronize();
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                               - t0)
                   .count();
      printf("rep %d: %.4f s (floats in, floats out)\n", rep, s);
      best = std::min(best, s);
    }
  double mpmac = double(rows) * chunks * N * (N + 1) / 2;
  printf("best %.4f s; %.3g MP-MAC -> %.3g MP-MAC/s\n", best, mpmac,
         mpmac / best);
  for(auto &[k, v] : MG.timings())
    printf("  phase %-22s %9.2f ms (all reps)\n", k.c_str(), v);
  for(auto &[k, v] : MG.gram().timings())
    printf("  phase %-22s %9.2f ms (all reps)\n", k.c_str(), v);

  // ---- check sampled entries against exact GMP sums
  std::vector<int64_t> cols;
  if(sample_cols == 0 || sample_cols >= N)
    for(int64_t j = 0; j < N; ++j)
      cols.push_back(j);
  else
    {
      std::set<int64_t> s = {0, N - 1, std::min<int64_t>(N - 1, 128),
                             std::min<int64_t>(N - 1, N / 2)};
      std::mt19937_64 rng(seed + 3);
      while((int)s.size() < sample_cols)
        s.insert(int64_t(rng() % N));
      cols.assign(s.begin(), s.end());
    }
  const size_t nc = cols.size();
  const int64_t K = reuse ? rows : rows * chunks;
  std::vector<__mpz_struct> A(nc * K);
  std::vector<int64_t> Ae(nc * K);
  std::vector<double> colmax(nc, -1e300);
  std::vector<uint32_t> w(64);
  for(size_t q = 0; q < nc; ++q)
    for(int64_t k = 0; k < K; ++k)
      {
        mpz_init(&A[q * K + k]);
        int64_t ex = gen_entry(seed, k / rows, k % rows, cols[q], L, w.data());
        Ae[q * K + k] = ex;
        if(ex != MPF_EZ)
          {
            tc_to_mpz(&A[q * K + k], w.data(), L);
            colmax[q] = std::max(colmax[q], log2_abs(&A[q * K + k], ex));
          }
      }
  std::vector<int32_t> hm(size_t(L) * total);
  std::vector<int64_t> he(total);
  cudaMemcpy(hm.data(), d_qm, hm.size() * 4, cudaMemcpyDeviceToHost);
  cudaMemcpy(he.data(), d_qe, he.size() * 8, cudaMemcpyDeviceToHost);
  mpz_t S, t, O;
  mpz_inits(S, t, O, nullptr);
  double worst_norm = -1e300, worst_rel = -1e300;
  long checked = 0, bad = 0;
  const double nbits = std::min(8.0 * D - 2, 32.0 * L - 2);
  for(size_t qj = 0; qj < nc; ++qj)
    for(size_t qi = qj; qi < nc; ++qi)
      {
        // exact sum S * 2^Emin
        int64_t Emin = INT64_MAX;
        for(int64_t k = 0; k < K; ++k)
          if(Ae[qi * K + k] != MPF_EZ && Ae[qj * K + k] != MPF_EZ)
            Emin = std::min(Emin, Ae[qi * K + k] + Ae[qj * K + k]);
        mpz_set_ui(S, 0);
        for(int64_t k = 0; k < K && Emin != INT64_MAX; ++k)
          if(Ae[qi * K + k] != MPF_EZ && Ae[qj * K + k] != MPF_EZ)
            {
              mpz_mul(t, &A[qi * K + k], &A[qj * K + k]);
              mpz_mul_2exp(t, t, Ae[qi * K + k] + Ae[qj * K + k] - Emin);
              mpz_add(S, S, t);
            }
        if(reuse)
          mpz_mul_ui(S, S, chunks);
        int64_t e = Gram_Accumulator::packed_index(N, cols[qi], cols[qj]);
        for(int q = 0; q < L; ++q)
          w[q] = uint32_t(hm[size_t(q) * total + e]);
        int64_t oe = he[e];
        if(oe == MPF_EZ)
          mpz_set_ui(O, 0);
        else
          tc_to_mpz(O, w.data(), L);
        // diff = O 2^oe - S 2^Emin, at the smaller exponent
        if(Emin == INT64_MAX)
          Emin = oe == MPF_EZ ? 0 : oe;
        int64_t ez = oe == MPF_EZ ? Emin : std::min(oe, Emin);
        mpz_mul_2exp(t, O, oe == MPF_EZ ? 0 : oe - ez);
        mpz_mul_2exp(S, S, Emin - ez);
        mpz_sub(t, t, S);
        double lerr = log2_abs(t, ez);
        double lq = log2_abs(S, ez);
        double lbound = colmax[qi] + colmax[qj] + std::log2(double(rows * chunks));
        double rn = lerr - lbound, rr = lerr - lq;
        worst_norm = std::max(worst_norm, rn);
        if(lq > -1e299)
          worst_rel = std::max(worst_rel, rr);
        ++checked;
        if(rn > -(nbits - 4))
          ++bad;
      }
  printf("CHECK %s: %ld entries; worst error / (rows max|P_i| max|P_j|) = 2^%.1f "
         "(gate 2^-%.0f), worst error / |Q_ij| = 2^%.1f\n",
         bad ? "FAILED" : "PASSED", checked, worst_norm, nbits - 4, worst_rel);
  return bad ? 1 : 0;
}
