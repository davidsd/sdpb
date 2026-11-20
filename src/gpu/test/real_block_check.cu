// GPU backend for SDPB: L^-1 B on real Schur blocks.  Reads the dim-5 Schur
// blocks S, the rows B and SDPviBe's own P = L^-1 B dumped by
// bench/sdpvibe_dump_schur.py (SDPviBe by Antoine Tilloy, used unchanged),
// pads the blocks to a multiple of 32 (identity), factors them with
// Block_Factor, forms P = X B, and reports
//   * per-column normwise difference from SDPviBe's P (all blocks, on the GPU),
//   * residuals |L L^T - S| (sampled, GMP),
//   * with --gmp-ref: block 0 against a GMP Cholesky and forward substitution at
//     2F+128 bits, for both the backend's and SDPviBe's P (sampled columns).
//   real_block_check --dir <dump> --nb 36 --P 945 --N 3406 [--gmp-ref --cols 4]

#include "../linalg/Block_Factor.hxx"

#include <gmp.h>

#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
  bool flag(int argc, char **argv, const char *name)
  {
    for(int i = 1; i < argc; ++i)
      if(!strcmp(argv[i], name))
        return true;
    return false;
  }
  // read SDPviBe's MPF (L, nb, R, C) row-major into a Mpf_Batch with rows = C,
  // cols = R (the transpose, column-major), batch nb
  void read_mpf(const std::string &base, Mpf_Batch &M)
  {
    const size_t ent = size_t(M.entries());
    std::vector<int32_t> hm(size_t(M.L) * ent);
    std::vector<int64_t> he(ent);
    FILE *f = fopen((base + ".m.bin").c_str(), "rb");
    if(!f || fread(hm.data(), 4, hm.size(), f) != hm.size())
      throw std::runtime_error("read " + base + ".m.bin");
    fclose(f);
    f = fopen((base + ".e.bin").c_str(), "rb");
    if(!f || fread(he.data(), 8, he.size(), f) != he.size())
      throw std::runtime_error("read " + base + ".e.bin");
    fclose(f);
    cudaMemcpy(M.m, hm.data(), hm.size() * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(M.e, he.data(), he.size() * 8, cudaMemcpyHostToDevice);
  }
  __global__ void fill_identity(int32_t *m, int64_t *e, int64_t total, int L,
                                int64_t n)
  {
    int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    if(idx >= total)
      return;
    int64_t rem = idx % (n * n), i = rem % n, j = rem / n;
    for(int t = 0; t < L; ++t)
      m[(int64_t)t * total + idx] = (i == j && t == L - 1) ? 0x20000000 : 0;
    e[idx] = (i == j) ? -(32LL * L - 3) : MPF_EZ;
  }
  // per (block, column): max over rows of e + bitlen(m) (nonzero entries)
  __global__ void col_top(const int32_t *m, const int64_t *e, int64_t total,
                          int L, int64_t R, int64_t C, long long *out)
  {
    int64_t q = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    int64_t batch = total / (R * C);
    if(q >= batch * C)
      return;
    int64_t b = q / C, j = q % C;
    long long best = LLONG_MIN;
    for(int64_t i = 0; i < R; ++i)
      {
        int64_t idx = b * R * C + j * R + i;
        if(e[idx] == MPF_EZ)
          continue;
        uint32_t ext = m[(int64_t)(L - 1) * total + idx] < 0 ? 0xffffffffu : 0u;
        int bl = 0;
        for(int t = L - 1; t >= 0; --t)
          {
            uint32_t w = uint32_t(m[(int64_t)t * total + idx]) ^ ext;
            if(w)
              {
                bl = 32 * t + 32 - __clz(w);
                break;
              }
          }
        best = max(best, (long long)(e[idx] + bl));
      }
    out[q] = best;
  }
  std::vector<long long> tops(const Mpf_Batch &A)
  {
    long long *d;
    const int64_t nq = int64_t(A.batch) * A.cols;
    cudaMalloc(&d, nq * 8);
    col_top<<<(nq + 255) / 256, 256>>>(A.m, A.e, A.entries(), A.L, A.rows,
                                       A.cols, d);
    std::vector<long long> h(nq);
    cudaMemcpy(h.data(), d, nq * 8, cudaMemcpyDeviceToHost);
    cudaFree(d);
    return h;
  }
  // one entry of a Mpf_Batch as mpf (downloads the entry's limbs)
  void entry(mpf_t out, const Mpf_Batch &A, int64_t idx)
  {
    std::vector<int32_t> w(A.L);
    int64_t ex;
    for(int t = 0; t < A.L; ++t)
      cudaMemcpy(&w[t], A.m + int64_t(t) * A.entries() + idx, 4,
                 cudaMemcpyDeviceToHost);
    cudaMemcpy(&ex, A.e + idx, 8, cudaMemcpyDeviceToHost);
    if(ex == MPF_EZ)
      {
        mpf_set_ui(out, 0);
        return;
      }
    mpz_t z;
    mpz_init(z);
    mpz_import(z, A.L, -1, 4, 0, 0, w.data());
    if(uint32_t(w[A.L - 1]) >> 31)
      {
        mpz_t p;
        mpz_init(p);
        mpz_ui_pow_ui(p, 2, 32 * A.L);
        mpz_sub(z, z, p);
        mpz_clear(p);
      }
    mpf_set_z(out, z);
    if(ex >= 0)
      mpf_mul_2exp(out, out, ex);
    else
      mpf_div_2exp(out, out, -ex);
    mpz_clear(z);
  }
  // block 0 of A (rows x cols) to host mpf, column-major
  std::vector<__mpf_struct> block0(const Mpf_Batch &A)
  {
    const int64_t n = A.rows * A.cols;
    std::vector<int32_t> hm(size_t(A.L) * n);
    std::vector<int64_t> he(n);
    for(int t = 0; t < A.L; ++t)
      cudaMemcpy(hm.data() + size_t(t) * n, A.m + int64_t(t) * A.entries(),
                 n * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(he.data(), A.e, n * 8, cudaMemcpyDeviceToHost);
    std::vector<__mpf_struct> out(n);
    std::vector<uint32_t> w(A.L);
    mpz_t z, p;
    mpz_inits(z, p, nullptr);
    mpz_ui_pow_ui(p, 2, 32 * A.L);
    for(int64_t q = 0; q < n; ++q)
      {
        mpf_init(&out[q]);
        if(he[q] == MPF_EZ)
          continue;
        for(int t = 0; t < A.L; ++t)
          w[t] = uint32_t(hm[size_t(t) * n + q]);
        mpz_import(z, A.L, -1, 4, 0, 0, w.data());
        if(w[A.L - 1] >> 31)
          mpz_sub(z, z, p);
        mpf_set_z(&out[q], z);
        if(he[q] >= 0)
          mpf_mul_2exp(&out[q], &out[q], he[q]);
        else
          mpf_div_2exp(&out[q], &out[q], -he[q]);
      }
    mpz_clears(z, p, nullptr);
    return out;
  }
  double lg(const mpf_t x)
  {
    if(mpf_sgn(x) == 0)
      return -1e9;
    long ex;
    double d = mpf_get_d_2exp(&ex, x);
    return std::log2(std::fabs(d)) + ex;
  }
}

int main(int argc, char **argv)
{
  const std::string dir = arg(argc, argv, "--dir", ".");
  const int nb = atoi(arg(argc, argv, "--nb", "36"));
  const int P = atoi(arg(argc, argv, "--P", "945"));
  const int N = atoi(arg(argc, argv, "--N", "3406"));
  const int L = atoi(arg(argc, argv, "--L", "28"));
  const int reps = atoi(arg(argc, argv, "--reps", "2"));
  const int ncols = atoi(arg(argc, argv, "--cols", "4"));
  // --chunk c: the dump's c-th chunk of blocks (S<c>, B<c>, P<c>)
  const std::string ch = arg(argc, argv, "--chunk", "0");
  const bool gref = flag(argc, argv, "--gmp-ref");
  Factor_Options fo;
  fo.equilibrate = atoi(arg(argc, argv, "--equilibrate", "1"));
  fo.panel_refine = atoi(arg(argc, argv, "--panel-refine", "1"));
  const double agreement_limit = atof(arg(argc, argv, "--agreement-log2", "-700"));
  const double residual_limit = atof(arg(argc, argv, "--residual-log2", "0"));
  const int n = (P + 31) / 32 * 32, D = 112;
  Gram_Options opt;
  opt.prime_group = 32;
  opt.fused_residues = true;
  for(int i = 1; i < argc; ++i)
    {
      if(!strcmp(argv[i], "--old-leaf"))
        opt.fused_leaf = false;
      if(!strcmp(argv[i], "--no-lookahead"))
        opt.leaf_lookahead = false;
      if(!strcmp(argv[i], "--crt256"))
        opt.crt_base128 = false;
      if(!strcmp(argv[i], "--no-reduce-crt"))
        opt.reduce_to_crt = false;
      if(!strcmp(argv[i], "--slow-rounding"))
        opt.fast_rounding = false;
    }
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("device %s; real Schur blocks from %s: %d blocks of %d (padded to %d), N=%d\n",
         prop.name, dir.c_str(), nb, P, n, N);
  Block_Factor F(n, nb, L, D, opt);
  // S (symmetric: SDPviBe's row-major block is our column-major block)
  Mpf_Batch S0 = F.alloc(P, P), Bt = F.alloc(N, P), Pt = F.alloc(N, P);
  read_mpf(dir + "/S" + ch, S0);
  read_mpf(dir + "/B" + ch, Bt);
  read_mpf(dir + "/P" + ch, Pt);
  Mpf_Batch B0 = F.alloc(P, N), Pref = F.alloc(P, N);
  mpf_transpose(Bt, B0, 0);
  mpf_transpose(Pt, Pref, 0);
  F.release(Bt);
  F.release(Pt);
  Mpf_Batch S = F.alloc(n, n), B = F.alloc(n, N);
  fill_identity<<<(S.entries() + 255) / 256, 256>>>(S.m, S.e, S.entries(), L, n);
  mpf_sub_block(S, 0, 0, S0, true, 0);
  mpf_zero(B, 0);
  mpf_sub_block(B, 0, 0, B0, true, 0);
  F.release(B0); // B's top P rows are B0
  cudaError_t ce = cudaDeviceSynchronize();
  if(ce != cudaSuccess)
    {
      printf("CUDA error %s\n", cudaGetErrorString(ce));
      return 2;
    }
  Mpf_Batch Lf = F.alloc(n, n), X = F.alloc(n, n), Pn;
  double bf = 1e30, bp = 1e30;
  for(int rep = 0; rep < reps; ++rep)
    {
      cudaDeviceSynchronize();
      auto t0 = std::chrono::steady_clock::now();
      F.factor(S, Lf, X, fo);
      cudaDeviceSynchronize();
      auto t1 = std::chrono::steady_clock::now();
      if(Pn.m)
        F.release(Pn);
      F.solve(Lf, X, B, Pn, 0);
      cudaDeviceSynchronize();
      auto t2 = std::chrono::steady_clock::now();
      double f = std::chrono::duration<double>(t1 - t0).count(),
             p = std::chrono::duration<double>(t2 - t1).count();
      printf("rep %d: factor %.4f s, P = X B %.4f s\n", rep, f, p);
      bf = std::min(bf, f);
      bp = std::min(bp, p);
    }
  printf("best: factor %.4f s + X B %.4f s = %.4f s\n", bf, bp, bf + bp);
  ce = cudaDeviceSynchronize();
  if(ce != cudaSuccess)
    {
      printf("CUDA error %s\n", cudaGetErrorString(ce));
      return 2;
    }
  // the backend's P (top P rows) against SDPviBe's
  F.clear_cache();
  F.release(X);
  Mpf_Batch Pb = F.alloc(P, N);
  mpf_sub_block(Pn, 0, 0, Pb, false, 0);
  F.release(Pn);
  Mpf_Batch Dm = F.alloc(P, N);
  mpf_batch_add(Pb, Pref, -1, Dm, 0);
  auto tD = tops(Dm), tR = tops(Pref);
  F.release(Dm);
  double worst = -1e9, sum = 0;
  long cnt = 0;
  for(size_t q = 0; q < tD.size(); ++q)
    {
      if(tR[q] == LLONG_MIN)
        continue;
      double r = (tD[q] == LLONG_MIN) ? -1e9 : double(tD[q] - tR[q]);
      worst = std::max(worst, r);
      if(r > -1e8)
        {
          sum += r;
          ++cnt;
        }
    }
  printf("backend P vs SDPviBe P, per column of every block: max_i|dP| / max_i|P| <= 2^%.0f "
         "(mean exponent %.1f over %ld columns)\n",
         worst, cnt ? sum / cnt : 0.0, cnt);
  // residuals |L L^T - S| (sampled, block 0..nb-1)
  mpf_set_default_prec(2 * 32 * L + 256);
  mpf_t a, b2, s, t;
  mpf_inits(a, b2, s, t, nullptr);
  double res = -1e9;
  srand(3);
  for(int q = 0; q < 12; ++q)
    {
      int64_t bb = rand() % nb, i = rand() % P, j = rand() % P;
      if(i < j)
        std::swap(i, j);
      mpf_set_ui(s, 0);
      for(int64_t k = 0; k <= j; ++k)
        {
          entry(a, Lf, bb * n * n + k * n + i);
          entry(b2, Lf, bb * n * n + k * n + j);
          mpf_mul(t, a, b2);
          mpf_add(s, s, t);
        }
      entry(a, S, bb * n * n + j * n + i);
      mpf_sub(s, s, a);
      entry(a, S, bb * n * n + i * n + i);
      entry(b2, S, bb * n * n + j * n + j);
      mpf_mul(a, a, b2);
      mpf_sqrt(a, a);
      res = std::max(res, lg(s) - lg(a));
    }
  printf("residual (12 sampled entries): max|LL^T - S|/sqrt(S_ii S_jj) = 2^%.1f\n", res);
  bool ok = worst < agreement_limit && res < residual_limit;
  // Agreement depends on conditioning and the comparator's own error. Keep the
  // historical default, but allow explicit gates for real Q (see test recipe).
  if(gref)
    {
      auto t0 = std::chrono::steady_clock::now();
      std::vector<__mpf_struct> R = block0(S0); // P x P, column-major
      // diagnostics: symmetry and diagonal range of block 0
      {
        double asym = -1e9, dmin = 1e9, dmax = -1e9;
        for(int64_t j = 0; j < P; j += 37)
          for(int64_t i = j + 1; i < P; i += 41)
            {
              mpf_sub(t, &R[j * P + i], &R[i * P + j]);
              asym = std::max(asym, lg(t) - lg(&R[j * P + i]));
            }
        for(int64_t i = 0; i < P; ++i)
          {
            dmin = std::min(dmin, lg(&R[i * P + i]));
            dmax = std::max(dmax, lg(&R[i * P + i]));
          }
        printf("block 0: relative asymmetry <= 2^%.1f, diagonal 2^%.1f .. 2^%.1f\n", asym, dmin, dmax);
      }
      for(int64_t j = 0; j < P; ++j)
        {
          mpf_set(s, &R[j * P + j]);
          for(int64_t k = 0; k < j; ++k)
            {
              mpf_mul(t, &R[k * P + j], &R[k * P + j]);
              mpf_sub(s, s, t);
            }
          if(mpf_sgn(s) <= 0)
            {
              printf("GMP Cholesky: pivot %ld not positive (2^%.1f)\n", long(j), lg(s));
              mpf_abs(s, s);
            }
          mpf_sqrt(&R[j * P + j], s);
          for(int64_t i = j + 1; i < P; ++i)
            {
              mpf_set(s, &R[j * P + i]);
              for(int64_t k = 0; k < j; ++k)
                {
                  mpf_mul(t, &R[k * P + i], &R[k * P + j]);
                  mpf_sub(s, s, t);
                }
              mpf_div(&R[j * P + i], s, &R[j * P + j]);
            }
        }
      double tch = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      double eb = -1e9, es = -1e9;
      std::vector<__mpf_struct> y(P);
      for(auto &v : y)
        mpf_init(&v);
      int used = 0;
      for(int c = 0; c < 4 * ncols && used < ncols; ++c)
        {
          const int64_t col = (int64_t(c) * 977) % N;
          double ymax = -1e9;
          for(int64_t i = 0; i < P; ++i)
            {
              entry(s, B, col * n + i);
              for(int64_t k = 0; k < i; ++k)
                {
                  mpf_mul(t, &R[k * P + i], &y[k]);
                  mpf_sub(s, s, t);
                }
              mpf_div(&y[i], s, &R[i * P + i]);
              ymax = std::max(ymax, lg(&y[i]));
            }
          if(c == 0)
            for(int64_t i = 0; i < 3; ++i)
              {
                entry(a, Pb, col * P + i);
                entry(b2, Pref, col * P + i);
                entry(s, B, col * n + i);
                printf("  col %ld row %ld: B 2^%.2f sign %d, y 2^%.2f sign %d, backend P 2^%.2f sign %d, "
                       "SDPviBe P 2^%.2f sign %d\n", long(col), long(i), lg(s), mpf_sgn(s), lg(&y[i]),
                       mpf_sgn(&y[i]), lg(a), mpf_sgn(a), lg(b2), mpf_sgn(b2));
              }
          if(ymax < -1e8)
            continue; // zero column of B: nothing to compare
          ++used;
          double cb = -1e9, cs = -1e9;
          int64_t wr = 0;
          for(int64_t i = 0; i < P; ++i)
            {
              entry(a, Pb, col * P + i);
              mpf_sub(a, a, &y[i]);
              if(lg(a) - ymax > cb)
                {
                  cb = lg(a) - ymax;
                  wr = i;
                }
              eb = std::max(eb, lg(a) - ymax);
              entry(a, Pref, col * P + i);
              mpf_sub(a, a, &y[i]);
              es = std::max(es, lg(a) - ymax);
              cs = std::max(cs, lg(a) - ymax);
            }
          entry(a, Pb, col * P + wr);
          entry(b2, Pref, col * P + wr);
          printf("  col %ld: max|y| 2^%.1f; worst row %ld: y 2^%.2f (sign %d), backend 2^%.2f (sign %d), "
                 "SDPviBe 2^%.2f (sign %d); column error backend 2^%.1f, SDPviBe 2^%.1f\n", long(col),
                 ymax, long(wr), lg(&y[wr]), mpf_sgn(&y[wr]), lg(a), mpf_sgn(a), lg(b2), mpf_sgn(b2),
                 cb, cs);
        }
      printf("block 0 vs GMP (Cholesky %.0f s, %d nonzero columns): backend max|P - P_exact|/max|P_exact| = "
             "2^%.1f, SDPviBe 2^%.1f\n", tch, used, eb, es);
    }
  printf("CHECK %s\n", ok ? "PASSED" : "FAILED");
  return ok ? 0 : 1;
}
