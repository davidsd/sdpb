// q_check: an independent check of the backend (Block_Factor) on real Schur
// blocks dumped by SDPviBe, with factor options, full factor/inverse export
// and timing controls. For every block it writes K selected columns of the
// backend's P = L^-1 B together with B and SDPviBe's P for the same columns, read
// straight from the dump files on the host (for ref_solve.c, an MPFR
// reference). Also: every column of every block, backend P against SDPviBe's P,
// computed on the host with GMP integers (not the backend's own comparison
// kernels). Selected columns per block: the W columns where the two codes
// disagree most, the backend's own columns (c*977 mod N, first 3 nonzero), the rest
// random nonzero.
//   q_check --dir <dump> --chunk 0 --nb 36 --P 945 --N 3406 --K 48 --W 16
//   --reps 4 --out cols.bin
#include "../linalg/Block_Factor.hxx"

#include <gmp.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sdpb_gpu;

namespace
{
  const char *arg(int argc, char **argv, const char *name, const char *dflt)
  {
    for(int i = 1; i + 1 < argc; ++i)
      if(!strcmp(argv[i], name))
        return argv[i + 1];
    return dflt;
  }
  template <class T>
  std::vector<T> slurp(const std::string &path, size_t count)
  {
    std::vector<T> v(count);
    FILE *f = fopen(path.c_str(), "rb");
    if(!f || fread(v.data(), sizeof(T), count, f) != count)
      throw std::runtime_error("read " + path);
    if(fgetc(f) != EOF)
      throw std::runtime_error("trailing data in " + path);
    fclose(f);
    return v;
  }
  const int64_t EZ = -(1LL << 40);
  // a host float: limb t at m[t * plane + idx], exponent e[idx]
  struct HostMpf
  {
    const int32_t *m;
    const int64_t *e;
    int64_t plane;
  };
  void to_mpz(mpz_t z, const HostMpf &A, int64_t idx, int L, mpz_t two)
  {
    uint32_t w[64];
    for(int t = 0; t < L; ++t)
      w[t] = uint32_t(A.m[t * A.plane + idx]);
    mpz_import(z, L, -1, 4, 0, 0, w);
    if(w[L - 1] >> 31)
      mpz_sub(z, z, two);
  }
} // namespace

int main(int argc, char **argv)
{
  const std::string dir = arg(argc, argv, "--dir", ".");
  const std::string ch = arg(argc, argv, "--chunk", "0");
  const int nb = atoi(arg(argc, argv, "--nb", "36"));
  const int P = atoi(arg(argc, argv, "--P", "945"));
  const int N = atoi(arg(argc, argv, "--N", "3406"));
  const int L = atoi(arg(argc, argv, "--L", "28"));
  const int K = atoi(arg(argc, argv, "--K", "48"));
  const int W = atoi(arg(argc, argv, "--W", "16"));
  const int reps = atoi(arg(argc, argv, "--reps", "4"));
  const std::string out = arg(argc, argv, "--out", "cols.bin");
  const int n = (P + 31) / 32 * 32,
            D = atoi(arg(argc, argv, "--digits", std::to_string(4 * L).c_str()));
  for(const char *flag : {"--panel-leaf", "--update-tile", "--auxiliary-digits"})
    if(atoi(arg(argc, argv, flag, "0")))
      throw std::invalid_argument(std::string(flag)
                                  + ": this experimental option is not included");
  Gram_Options opt; // the defaults of real_block_check
  opt.prime_group = 32;
  opt.fused_residues = true;
  opt.crt_guard_bits = atoi(arg(argc, argv, "--crt-guard", "0"));
  opt.share_workspace = atoi(arg(argc, argv, "--share-workspace", "1"));
  auto T0 = std::chrono::steady_clock::now();
  const int64_t sp = int64_t(nb) * P * P, bp = int64_t(nb) * P * N;
  auto Sm = slurp<int32_t>(dir + "/S" + ch + ".m.bin", size_t(L) * sp);
  auto Se = slurp<int64_t>(dir + "/S" + ch + ".e.bin", sp);
  auto Bm = slurp<int32_t>(dir + "/B" + ch + ".m.bin", size_t(L) * bp);
  auto Be = slurp<int64_t>(dir + "/B" + ch + ".e.bin", bp);
  auto Pm = slurp<int32_t>(dir + "/P" + ch + ".m.bin", size_t(L) * bp);
  auto Pe = slurp<int64_t>(dir + "/P" + ch + ".e.bin", bp);
  printf("read %s chunk %s (nb=%d P=%d N=%d L=%d) in %.1f s\n", dir.c_str(),
         ch.c_str(), nb, P, N, L,
         std::chrono::duration<double>(std::chrono::steady_clock::now() - T0)
           .count());
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("device %s\n", prop.name);

  if(atoi(arg(argc, argv, "--retain-pool", "0")))
    {
      cudaMemPool_t pool;
      cudaDeviceGetDefaultMemPool(&pool, 0);
      uint64_t bytes = 64ULL << 30;
      cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &bytes);
    }
  Block_Factor F(n, nb, L, D, opt);
  // S: the dump's row-major P x P blocks, uploaded as column-major (= S^T = S;
  // ref_solve checks that the stored blocks are exactly symmetric)
  Mpf_Batch S0 = F.alloc(P, P);
  cudaMemcpy(S0.m, Sm.data(), Sm.size() * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(S0.e, Se.data(), Se.size() * 8, cudaMemcpyHostToDevice);
  // B: the dump's row-major P x N is a column-major N x P (= B^T); transpose
  Mpf_Batch Bt = F.alloc(N, P);
  cudaMemcpy(Bt.m, Bm.data(), Bm.size() * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(Bt.e, Be.data(), Be.size() * 8, cudaMemcpyHostToDevice);
  Mpf_Batch B0 = F.alloc(P, N);
  mpf_transpose(Bt, B0, 0);
  F.release(Bt);
  // pad to n with the identity (S) and zero rows (B)
  Mpf_Batch S = F.alloc(n, n), B = F.alloc(n, N);
  {
    std::vector<int32_t> hm(size_t(L) * S.entries(), 0);
    std::vector<int64_t> he(S.entries(), EZ);
    for(int b = 0; b < nb; ++b)
      for(int i = 0; i < n; ++i)
        {
          const int64_t idx = int64_t(b) * n * n + int64_t(i) * n + i;
          hm[size_t(L - 1) * S.entries() + idx]
            = 1 << 29; // 2^(32L-3) * 2^-(32L-3) = 1
          he[idx] = -(32LL * L - 3);
        }
    cudaMemcpy(S.m, hm.data(), hm.size() * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(S.e, he.data(), he.size() * 8, cudaMemcpyHostToDevice);
  }
  mpf_sub_block(S, 0, 0, S0, true, 0);
  mpf_zero(B, 0);
  mpf_sub_block(B, 0, 0, B0, true, 0);
  F.release(S0);
  F.release(B0);
  if(cudaDeviceSynchronize() != cudaSuccess)
    return printf("CUDA error after upload\n"), 2;

  Mpf_Batch Lf = F.alloc(n, n), X = F.alloc(n, n), Pn;
  std::vector<double> tf, tp;
  for(int rep = 0; rep < reps; ++rep)
    {
      cudaDeviceSynchronize();
      auto t0 = std::chrono::steady_clock::now();
      Factor_Options fo{atoi(arg(argc, argv, "--equilibrate", "1")) != 0,
                        atoi(arg(argc, argv, "--panel-refine", "1")),
                        atoi(arg(argc, argv, "--refine-min", "0")),
                        atoi(arg(argc, argv, "--correction-digits", "0")),
                        atoi(arg(argc, argv, "--gram-min", "0")),
                        atoi(arg(argc, argv, "--skip-triangular-zeros", "1")) != 0};
      fo.lower_update = atoi(arg(argc, argv, "--lower-update", "1"));
      fo.correction_limit = atoi(arg(argc, argv, "--correction-limit", "0"));
      fo.structured_solve = atoi(arg(argc, argv, "--structured-solve", "1"));
      fo.reuse_residues = atoi(arg(argc, argv, "--reuse-residues", "1"));
      fo.refine_max_n = atoi(arg(argc, argv, "--refine-max", "0"));
      fo.corrupt_leaf_limbs = atoi(arg(argc, argv, "--corrupt-leaf-limbs", "0"));
      fo.log_panels = rep == 0 && atoi(arg(argc, argv, "--log-panels", "0"));
      F.factor(S, Lf, X, fo);
      printf("POLICY rep %d: panels uncorrected %d, restarts %d\n", rep,
             F.panels_skipped(), F.correction_restarts());
      cudaDeviceSynchronize();
      auto t1 = std::chrono::steady_clock::now();
      if(Pn.m)
        F.release(Pn);
      F.solve(Lf, X, B, Pn, atoi(arg(argc, argv, "--solve-refine", "0")));
      cudaDeviceSynchronize();
      auto t2 = std::chrono::steady_clock::now();
      tf.push_back(std::chrono::duration<double>(t1 - t0).count());
      tp.push_back(std::chrono::duration<double>(t2 - t1).count());
      printf("TIME rep %d: factor %.4f s, P = X B %.4f s, total %.4f s\n", rep,
             tf.back(), tp.back(), tf.back() + tp.back());
    }
  F.check();
  if(cudaDeviceSynchronize() != cudaSuccess)
    return printf("CUDA error after solve\n"), 2;
  // the backend's P (n x N per block, column-major), transposed on the device
  // to N x n column-major = the dump's row-major layout (row i at i*N), to the
  // host
  printf("backend P: %ld x %ld, batch %d, L %d\n", long(Pn.rows),
         long(Pn.cols), Pn.batch, Pn.L);
  if(Pn.rows != n || Pn.cols != N || Pn.L != L || Pn.batch != nb)
    return printf("unexpected shape of P\n"), 2;
  F.clear_cache();
  auto export_matrix = [&](const char *flag, const Mpf_Batch &M) {
    const std::string factor_out = arg(argc, argv, flag, "");
    if(factor_out.empty())
      return;
    std::vector<int32_t> lm(size_t(L) * M.entries());
    std::vector<int64_t> le(M.entries());
    cudaMemcpy(lm.data(), M.m, lm.size() * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(le.data(), M.e, le.size() * 8, cudaMemcpyDeviceToHost);
    FILE *fm = fopen((factor_out + ".m.bin").c_str(), "wb"),
         *fe = fopen((factor_out + ".e.bin").c_str(), "wb");
    if(!fm || !fe)
      throw std::runtime_error("factor output");
    fwrite(lm.data(), 4, lm.size(), fm);
    fwrite(le.data(), 8, le.size(), fe);
    fclose(fm);
    fclose(fe);
  };
  export_matrix("--factor-out", Lf);
  export_matrix("--inverse-out", X);
  F.release(X);
  F.release(Lf);
  Mpf_Batch Pt = F.alloc(N, n);
  mpf_transpose(Pn, Pt, 0);
  if(cudaDeviceSynchronize() != cudaSuccess)
    return printf("CUDA error in transpose\n"), 2;
  std::vector<int32_t> Qm(size_t(L) * Pt.entries());
  std::vector<int64_t> Qe(Pt.entries());
  cudaMemcpy(Qm.data(), Pt.m, Qm.size() * 4, cudaMemcpyDeviceToHost);
  cudaMemcpy(Qe.data(), Pt.e, Qe.size() * 8, cudaMemcpyDeviceToHost);
  F.release(Pt);
  const HostMpf Hb{Qm.data(), Qe.data(), int64_t(Qe.size())},
    Hs{Pm.data(), Pe.data(), bp}, HB{Bm.data(), Be.data(), bp};
  // index of entry (row i, column j) of block b, both row-major
  auto ib = [&](int b, int i, int j) {
    return int64_t(b) * n * N + int64_t(i) * N + j;
  };
  auto is = [&](int b, int i, int j) {
    return int64_t(b) * P * N + int64_t(i) * N + j;
  };
  // rows P..n-1 of the backend's P must be zero (zero rows of B under the
  // identity)
  long padnz = 0;
  for(int b = 0; b < nb; ++b)
    for(int i = P; i < n; ++i)
      for(int j = 0; j < N; ++j)
        padnz += Qe[ib(b, i, j)] != EZ;
  // every column: log2 max_i |Pb - Psv| - log2 max_i |Psv| (GMP integers, exact
  // differences; the log2 are bit lengths, so within 1 of the true ratio)
  std::vector<double> dis(size_t(nb) * N, -1e9);
  std::vector<char> nz(size_t(nb) * N, 0);
  auto T1 = std::chrono::steady_clock::now();
#pragma omp parallel
  {
    mpz_t za, zb, two;
    mpz_inits(za, zb, two, nullptr);
    mpz_setbit(two, 32 * L);
    std::vector<long long> dtop(N), ptop(N);
#pragma omp for schedule(dynamic, 1)
    for(int b = 0; b < nb; ++b)
      {
        std::fill(dtop.begin(), dtop.end(), LLONG_MIN);
        std::fill(ptop.begin(), ptop.end(), LLONG_MIN);
        for(int i = 0; i < P; ++i)
          for(int j = 0; j < N; ++j)
            {
              const int64_t xs = is(b, i, j), xb = ib(b, i, j);
              const int64_t es = Pe[xs], eb = Qe[xb];
              if(es == EZ && eb == EZ)
                continue;
              if(es != EZ)
                {
                  to_mpz(za, Hs, xs, L, two);
                  ptop[j] = std::max(ptop[j],
                                     (long long)(es + mpz_sizeinbase(za, 2)));
                }
              else
                mpz_set_ui(za, 0);
              if(eb != EZ)
                to_mpz(zb, Hb, xb, L, two);
              else
                mpz_set_ui(zb, 0);
              // align to the smaller exponent and subtract exactly
              const int64_t e0 = (es == EZ)   ? eb
                                 : (eb == EZ) ? es
                                              : std::min(es, eb);
              if(es != EZ && es > e0)
                mpz_mul_2exp(za, za, es - e0);
              if(eb != EZ && eb > e0)
                mpz_mul_2exp(zb, zb, eb - e0);
              mpz_sub(za, za, zb);
              if(mpz_sgn(za))
                dtop[j]
                  = std::max(dtop[j], (long long)(e0 + mpz_sizeinbase(za, 2)));
            }
        for(int j = 0; j < N; ++j)
          if(ptop[j] != LLONG_MIN)
            {
              nz[size_t(b) * N + j] = 1;
              dis[size_t(b) * N + j]
                = (dtop[j] == LLONG_MIN) ? -1e9 : double(dtop[j] - ptop[j]);
            }
      }
    mpz_clears(za, zb, two, nullptr);
  }
  double worst = -1e9, sum = 0;
  long cnt = 0, ncols = 0, exact = 0;
  std::vector<double> bw(nb, -1e9);
  for(int64_t q = 0; q < int64_t(nb) * N; ++q)
    if(nz[q])
      {
        ++ncols;
        bw[q / N] = std::max(bw[q / N], dis[q]);
        worst = std::max(worst, dis[q]);
        if(dis[q] > -1e8)
          sum += dis[q], ++cnt;
        else
          ++exact;
      }
  printf("AGREE all %ld nonzero columns: max_i|Pb - Psv| / max_i|Psv| <= "
         "2^%.0f (bit-length bound), mean 2^%.1f; "
         "%ld columns identical; %ld nonzero padding entries; host compare "
         "%.1f s\n",
         ncols, worst, cnt ? sum / cnt : 0.0, exact, padnz,
         std::chrono::duration<double>(std::chrono::steady_clock::now() - T1)
           .count());
  for(int b = 0; b < nb; ++b)
    printf("AGREE block %d: worst 2^%.0f\n", b, bw[b]);
  // column selection
  std::vector<int32_t> col(size_t(nb) * K, -1), kind(size_t(nb) * K, 0);
  std::mt19937_64 rng(321 + atoi(ch.c_str()));
  for(int b = 0; b < nb; ++b)
    {
      std::vector<int> nzc, chosen;
      for(int j = 0; j < N; ++j)
        if(nz[size_t(b) * N + j])
          nzc.push_back(j);
      std::vector<int> byd = nzc;
      std::sort(byd.begin(), byd.end(), [&](int x, int y) {
        return dis[size_t(b) * N + x] > dis[size_t(b) * N + y];
      });
      auto take = [&](int j, int k) {
        if(std::find(chosen.begin(), chosen.end(), j) != chosen.end()
           || int(chosen.size()) >= K)
          return;
        kind[size_t(b) * K + chosen.size()] = k;
        col[size_t(b) * K + chosen.size()] = j;
        chosen.push_back(j);
      };
      for(int w = 0; w < W && w < int(byd.size()); ++w)
        take(byd[w], 1);
      int used = 0;
      for(int c = 0; c < 12 && used < 3; ++c)
        {
          int j = int((int64_t(c) * 977) % N);
          if(nz[size_t(b) * N + j])
            take(j, 2), ++used;
        }
      std::shuffle(nzc.begin(), nzc.end(), rng);
      for(int j : nzc)
        take(j, 0);
    }
  FILE *f = fopen(out.c_str(), "wb");
  fwrite("QCOL", 1, 4, f);
  int32_t h[5] = {nb, P, N, L, K};
  fwrite(h, 4, 5, f);
  fwrite(col.data(), 4, col.size(), f);
  fwrite(kind.data(), 4, kind.size(), f);
  std::vector<char> buf(4 * size_t(L) + 8);
  auto put = [&](const HostMpf &A, int64_t idx) {
    for(int t = 0; t < L; ++t)
      memcpy(buf.data() + 4 * t, &A.m[t * A.plane + idx], 4);
    memcpy(buf.data() + 4 * L, &A.e[idx], 8);
    fwrite(buf.data(), 1, buf.size(), f);
  };
  for(int b = 0; b < nb; ++b)
    for(int k = 0; k < K; ++k)
      {
        const int j = std::max(0, col[size_t(b) * K + k]);
        for(int i = 0; i < P; ++i)
          put(HB, is(b, i, j));
        for(int i = 0; i < P; ++i)
          put(Hb, ib(b, i, j));
        for(int i = 0; i < P; ++i)
          put(Hs, is(b, i, j));
      }
  fclose(f);
  std::vector<double> s1 = tf, s2 = tp, s3;
  for(size_t r = 1; r < tf.size(); ++r)
    s3.push_back(tf[r] + tp[r]);
  std::sort(s3.begin(), s3.end());
  if(!s3.empty())
    printf("TIME after warm-up (%zu reps): total min %.4f, median %.4f, max "
           "%.4f s\n",
           s3.size(), s3.front(), s3[s3.size() / 2], s3.back());
  printf("wrote %s; done in %.1f s\n", out.c_str(),
         std::chrono::duration<double>(std::chrono::steady_clock::now() - T0)
           .count());
  printf("CACHE_ALLOCATION_EVICTIONS %zu\n", F.cache_evictions());
  return 0;
}
