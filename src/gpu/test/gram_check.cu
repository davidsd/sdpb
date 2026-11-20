// GPU backend for SDPB: exactness check and benchmark of Gram_Accumulator.
//
// Generates a big-integer matrix A (rows x N per chunk, balanced base-256
// digits) from a counter-based hash, on the device for the GPU and again on
// the host for the reference, accumulates Q = sum_chunks A^T A on the GPU,
// and compares sampled (or all) entries of the lower triangle with GMP
// (schoolbook mpz products, independent of the residue/CRT path).
//
//   gram_check --N 3406 --rows 17520 --chunks 4 --digits 108 --sample-cols 40
//              [--mode random|allneg|sparse] [--group 16] [--base 256]
//              [--reps 2] [--seed 1]
//
// Prints per-phase times (with --timing), the INT8 rate and the check result.

#include "../rns/Gram_Accumulator.hxx"

#include <gmp.h>

#include <chrono>
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

  // mode 0: uniform digits; 1: every digit -128 (largest |a|); 2: random,
  // with 3/4 of the entries zero and every 7th column zero
  __host__ __device__ inline int8_t gen_digit(uint64_t seed, int mode,
                                              int64_t chunk, int64_t col,
                                              int64_t row, int t)
  {
    if(mode == 1)
      return -128;
    uint64_t h = splitmix(seed ^ splitmix(uint64_t(chunk) * 0x100000001B3ULL
                                          ^ splitmix(uint64_t(col) << 20
                                                     ^ uint64_t(row))));
    if(mode == 2)
      {
        if(col % 7 == 3 || (h & 3) != 0)
          return 0;
      }
    return (int8_t)(splitmix(h + uint64_t(t)) & 255);
  }

  __global__ void fill(int8_t *d, int64_t rows, int64_t N, int NB,
                       uint64_t seed, int mode, int64_t chunk)
  {
    int64_t e = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    if(e >= rows * N)
      return;
    int64_t col = e / rows, row = e % rows;
    for(int t = 0; t < NB; ++t)
      d[e * NB + t] = gen_digit(seed, mode, chunk, col, row, t);
  }

  void entry_mpz(mpz_t a, uint64_t seed, int mode, int64_t chunk,
                 int64_t col, int64_t row, int NB)
  {
    mpz_set_si(a, gen_digit(seed, mode, chunk, col, row, NB - 1));
    for(int t = NB - 2; t >= 0; --t)
      {
        mpz_mul_2exp(a, a, 8);
        long d = gen_digit(seed, mode, chunk, col, row, t);
        if(d >= 0)
          mpz_add_ui(a, a, d);
        else
          mpz_sub_ui(a, a, -d);
      }
  }

  // two's complement limbs -> mpz
  void limbs_mpz(mpz_t x, const std::vector<uint32_t> &w)
  {
    mpz_import(x, w.size(), -1, 4, 0, 0, w.data());
    if(w.back() >> 31)
      {
        mpz_t t;
        mpz_init(t);
        mpz_ui_pow_ui(t, 2, 32 * w.size());
        mpz_sub(x, x, t);
        mpz_clear(t);
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
  const int64_t N = atoll(arg(argc, argv, "--N", "200"));
  const int64_t rows = atoll(arg(argc, argv, "--rows", "1024"));
  const int chunks = atoi(arg(argc, argv, "--chunks", "2"));
  const int NB = atoi(arg(argc, argv, "--digits", "108"));
  const std::string mode_s = arg(argc, argv, "--mode", "random");
  const int mode = mode_s == "allneg" ? 1 : mode_s == "sparse" ? 2 : 0;
  const int sample_cols = atoi(arg(argc, argv, "--sample-cols", "0"));
  const int reps = atoi(arg(argc, argv, "--reps", "1"));
  const uint64_t seed = strtoull(arg(argc, argv, "--seed", "1"), nullptr, 10);
  // --reuse: every chunk is the same matrix (chunk 0), as a benchmark that
  // needs one chunk buffer; the GMP reference uses chunk 0 for every chunk
  const bool reuse = flag(argc, argv, "--reuse");
  Gram_Options opt;
  opt.prime_group = atoi(arg(argc, argv, "--group", "16"));
  opt.base = atoi(arg(argc, argv, "--base", "256"));
  opt.max_rect = atoi(arg(argc, argv, "--max-rect", "2048"));
  opt.timing = flag(argc, argv, "--timing");
  opt.crt_batch = atoll(arg(argc, argv, "--crt-batch", "524288"));
  opt.fused_residues = flag(argc, argv, "--fused");
  opt.fused_crt = flag(argc, argv, "--fused-crt");

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("device %s, %.1f GB; N=%lld rows/chunk=%lld chunks=%d digits=%d "
         "(%d bits) mode=%s group=%d base=%d%s\n",
         prop.name, prop.totalGlobalMem / 1e9, (long long)N, (long long)rows,
         chunks, NB, 8 * NB - 1, mode_s.c_str(), opt.prime_group, opt.base,
         opt.fused_residues ? " fused-residues" : "");

  // chunk buffers: one per chunk if memory allows, else regenerate per chunk
  const size_t chunk_bytes = size_t(rows) * N * NB;
  std::vector<int8_t *> d_chunks(chunks, nullptr);
  size_t free_b = 0, total_b = 0;
  cudaMemGetInfo(&free_b, &total_b);
  const bool resident = reuse || chunk_bytes * chunks < free_b * 0.55;
  const int nbuf = reuse ? 1 : resident ? chunks : 1;
  for(int c = 0; c < nbuf; ++c)
    if(cudaMalloc(&d_chunks[c], chunk_bytes) != cudaSuccess)
      {
        fprintf(stderr, "cannot allocate %.1f GB for chunk %d\n",
                chunk_bytes / 1e9, c);
        return 2;
      }
  auto gen = [&](int c) {
    int8_t *d = d_chunks[resident ? c : 0];
    int64_t ent = rows * N;
    fill<<<(ent + 255) / 256, 256>>>(d, rows, N, NB, seed, mode, c);
    return d;
  };
  if(resident)
    for(int c = 0; c < nbuf; ++c)
      gen(c);
  cudaDeviceSynchronize();

  Gram_Accumulator G(N, NB, rows * chunks, opt);
  const int W = G.limbs_out();
  printf("primes %zu (padded %d), M %d bits, CRT digits %d, output limbs %d; "
         "device buffers %.2f GB; chunks %s\n",
         G.basis().primes.size(), G.basis().num_padded, G.basis().M_bits,
         G.basis().crt_digits, W, G.device_bytes() / 1e9,
         reuse ? "one chunk reused" : resident ? "resident"
                                    : "regenerated per chunk (timing includes fill)");
  int32_t *d_out = nullptr;
  if(cudaMalloc(&d_out, size_t(W) * G.packed_size() * 4) != cudaSuccess)
    {
      fprintf(stderr, "cannot allocate output\n");
      return 2;
    }

  double best = 1e30;
  for(int rep = 0; rep < reps; ++rep)
    {
      double ops0 = G.int8_ops();
      cudaDeviceSynchronize();
      auto t0 = std::chrono::steady_clock::now();
      for(int c = 0; c < chunks; ++c)
        G.add(reuse ? d_chunks[0] : resident ? d_chunks[c] : gen(c), rows);
      G.finish(d_out);
      cudaDeviceSynchronize();
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                               - t0)
                   .count();
      double ops = G.int8_ops() - ops0;
      printf("rep %d: %.4f s total, %.3g INT8 ops issued, %.0f TOPS effective\n",
             rep, s, ops, ops / s / 1e12);
      best = std::min(best, s);
    }
  double mpmac = double(rows) * chunks * N * (N + 1) / 2;
  printf("best %.4f s; %.3g multiprecision multiply-adds (lower triangle) -> "
         "%.3g MP-MAC/s\n",
         best, mpmac, mpmac / best);
  for(auto &[k, v] : G.timings())
    printf("  phase %-22s %9.2f ms (all reps)\n", k.c_str(), v);

  // ---- check against GMP
  std::vector<int64_t> cols;
  if(sample_cols == 0 || sample_cols >= N)
    for(int64_t j = 0; j < N; ++j)
      cols.push_back(j);
  else
    {
      std::set<int64_t> s = {0, N - 1, std::min<int64_t>(N - 1, 31),
                             std::min<int64_t>(N - 1, 32),
                             std::min<int64_t>(N - 1, 127),
                             std::min<int64_t>(N - 1, 128),
                             std::min<int64_t>(N - 1, 255),
                             std::min<int64_t>(N - 1, 256),
                             std::min<int64_t>(N - 1, N / 2),
                             std::min<int64_t>(N - 1, N / 2 + 1)};
      std::mt19937_64 rng(seed + 17);
      while((int)s.size() < sample_cols)
        s.insert(int64_t(rng() % N));
      cols.assign(s.begin(), s.end());
    }
  const size_t nc = cols.size();
  const int64_t K = reuse ? rows : rows * chunks; // distinct rows
  printf("checking %zu x %zu sampled entries (all pairs of %zu columns) "
         "against GMP...\n",
         nc * (nc + 1) / 2, nc * (nc + 1) / 2, nc);
  fflush(stdout);
  std::vector<__mpz_struct> A(nc * K);
  for(size_t q = 0; q < nc; ++q)
    for(int64_t k = 0; k < K; ++k)
      {
        mpz_init(&A[q * K + k]);
        entry_mpz(&A[q * K + k], seed, mode, k / rows, cols[q], k % rows, NB);
      }
  const int64_t P = G.packed_size();
  std::vector<int32_t> h_out(size_t(W) * P);
  cudaMemcpy(h_out.data(), d_out, h_out.size() * 4, cudaMemcpyDeviceToHost);
  mpz_t ref, got, tmp;
  mpz_inits(ref, got, tmp, nullptr);
  long bad = 0, checked = 0;
  size_t maxbits = 0;
  std::vector<uint32_t> w(W);
  for(size_t qj = 0; qj < nc; ++qj)
    for(size_t qi = qj; qi < nc; ++qi)
      {
        int64_t i = cols[qi], j = cols[qj];
        mpz_set_ui(ref, 0);
        for(int64_t k = 0; k < K; ++k)
          mpz_addmul(ref, &A[qi * K + k], &A[qj * K + k]);
        if(reuse)
          mpz_mul_ui(ref, ref, chunks);
        int64_t e = Gram_Accumulator::packed_index(N, i, j);
        for(int t = 0; t < W; ++t)
          w[t] = uint32_t(h_out[size_t(t) * P + e]);
        limbs_mpz(got, w);
        maxbits = std::max(maxbits, mpz_sizeinbase(ref, 2));
        ++checked;
        if(mpz_cmp(ref, got) != 0)
          {
            if(bad < 5)
              gmp_printf("MISMATCH (%lld,%lld): ref %Zd\n got %Zd\n",
                         (long long)i, (long long)j, ref, got);
            ++bad;
          }
      }
  printf("CHECK %s: %ld of %ld entries differ from GMP (largest |Q_ij| %zu "
         "bits; bound %d bits)\n",
         bad ? "FAILED" : "PASSED", bad, checked, maxbits,
         G.basis().M_bits - 2);
  return bad ? 1 : 0;
}
