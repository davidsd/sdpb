// GPU backend for SDPB: check and benchmark of Block_Factor (batched Cholesky
// with explicit inverses, P = L^-1 B as products, optional refinement).
//
// SPD blocks S_b = M_b M_b^T (+ n I), M_b random floats with column scales
// 2^(-c j) (condition number up to ~2^(2 c n)), formed on the GPU; B random.
// Checks, in exact GMP arithmetic from the floats the GPU used and produced:
//   residuals  max |S - L L^T| / max|S|,  max |L P - B| / (max|L| max|P|)  (any n)
//   forward    against a GMP mpf Cholesky and substitution at 2F+128 bits:
//              max |L - L_ref| / max|L_ref|, max |P - P_ref| / max|P_ref|
//              (when --ref, small n)
//   block_factor_check --n 128 --ncols 64 --batch 2 --scale 0.5 --ref
//   block_factor_check --n 960 --ncols 3406 --batch 36 --reps 2 --check-entries 50

#include "../linalg/Block_Factor.hxx"

#include <gmp.h>

#include <chrono>
#include <cmath>
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
  // random normalized float, magnitude ~ 2^(-cs), sign random
  __global__ void fill_rand(int32_t *m, int64_t *e, int64_t total, int L,
                            int64_t rows, int64_t cols, uint64_t seed,
                            double cscale, double rspan)
  {
    int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    if(idx >= total)
      return;
    int64_t j = (idx % (rows * cols)) / rows;
    uint64_t h = splitmix(seed ^ splitmix(uint64_t(idx)));
    uint32_t w[64];
    for(int t = 0; t < L - 1; ++t)
      w[t] = uint32_t(splitmix(h + 1 + t));
    w[L - 1] = (uint32_t(splitmix(h + 99)) & 0x1FFFFFFFu) | 0x20000000u;
    if(h & 1)
      {
        unsigned long long c = 1;
        for(int t = 0; t < L; ++t)
          {
            unsigned long long v = (unsigned long long)(~w[t]) + c;
            w[t] = uint32_t(v);
            c = v >> 32;
          }
      }
    for(int t = 0; t < L; ++t)
      m[(int64_t)t * total + idx] = int32_t(w[t]);
    e[idx] = -(32LL * L - 2) - (long long)llround(cscale * double(j))
             - (long long)llround(rspan * double(idx % rows) / (rows > 1 ? rows - 1 : 1));
  }
  // identity * 2^s on the diagonal (adds n I to make S comfortably SPD)
  __global__ void fill_diag(int32_t *m, int64_t *e, int64_t total, int L,
                            int64_t n, int64_t s)
  {
    int64_t idx = blockIdx.x * (int64_t)blockDim.x + threadIdx.x;
    if(idx >= total)
      return;
    int64_t rem = idx % (n * n), i = rem % n, j = rem / n;
    for(int t = 0; t < L; ++t)
      m[(int64_t)t * total + idx] = (i == j && t == L - 1) ? 0x20000000 : 0;
    e[idx] = (i == j) ? s - (32LL * L - 3) : MPF_EZ;
  }

  void to_mpz(mpz_t x, const std::vector<int32_t> &hm, int64_t total, int L,
              int64_t idx)
  {
    std::vector<uint32_t> w(L);
    for(int t = 0; t < L; ++t)
      w[t] = uint32_t(hm[size_t(t) * total + idx]);
    mpz_import(x, L, -1, 4, 0, 0, w.data());
    if(w[L - 1] >> 31)
      {
        mpz_t p;
        mpz_init(p);
        mpz_ui_pow_ui(p, 2, 32 * L);
        mpz_sub(x, x, p);
        mpz_clear(p);
      }
  }
  struct Host
  {
    std::vector<int32_t> m;
    std::vector<int64_t> e;
    int64_t total;
    int L;
  };
  Host down(const Mpf_Batch &A)
  {
    Host H;
    H.total = A.entries();
    H.L = A.L;
    H.m.resize(size_t(A.L) * H.total);
    H.e.resize(H.total);
    cudaMemcpy(H.m.data(), A.m, H.m.size() * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(H.e.data(), A.e, H.e.size() * 8, cudaMemcpyDeviceToHost);
    return H;
  }
  void val(mpf_t out, const Host &H, int64_t idx)
  {
    if(H.e[idx] == MPF_EZ)
      {
        mpf_set_ui(out, 0);
        return;
      }
    mpz_t z;
    mpz_init(z);
    to_mpz(z, H.m, H.total, H.L, idx);
    mpf_set_z(out, z);
    if(H.e[idx] >= 0)
      mpf_mul_2exp(out, out, H.e[idx]);
    else
      mpf_div_2exp(out, out, -H.e[idx]);
    mpz_clear(z);
  }
  double lg(const mpf_t x)
  {
    if(mpf_sgn(x) == 0)
      return -1e9;
    long ex;
    double d = mpf_get_d_2exp(&ex, x);
    return std::log2(std::fabs(d)) + ex;
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
  const int n = atoi(arg(argc, argv, "--n", "128"));
  const int nc = atoi(arg(argc, argv, "--ncols", "64"));
  const int batch = atoi(arg(argc, argv, "--batch", "2"));
  const double cs = atof(arg(argc, argv, "--scale", "0.5"));
  const int reps = atoi(arg(argc, argv, "--reps", "1"));
  const int refine = atoi(arg(argc, argv, "--refine", "1"));
  const int checks = atoi(arg(argc, argv, "--check-entries", "0"));
  const bool ref = flag(argc, argv, "--ref");
  // --solve-q (with --ref): also Z = X^T (X B) = S^-1 B, the solve with S through
  // the explicit inverse, against GMP forward and backward substitution
  const bool solve_q = flag(argc, argv, "--solve-q");
  // --drop-cache: free the factorization's product objects before the solve
  // (memory at large batch; the next factor() re-creates them)
  const bool drop = flag(argc, argv, "--drop-cache");
  // --perturb-bits p: negative control. The GPU factors S + 2^-p S_00 E (E: ones
  // on the diagonal), while the checks use the unperturbed S; they must fail.
  const int perturb = atoi(arg(argc, argv, "--perturb-bits", "0"));
  const double row_span = atof(arg(argc, argv, "--row-span", "0"));
  const double forward_limit = atof(arg(argc, argv, "--forward-log2", "0"));
  const int ridge = atoi(arg(argc, argv, "--ridge-bits", "40"));
  Factor_Options fo;
  fo.equilibrate = atoi(arg(argc, argv, "--equilibrate", "1"));
  fo.panel_refine = atoi(arg(argc, argv, "--panel-refine", "1"));
  const int L = atoi(arg(argc, argv, "--L", "28")), D = 4 * L;
  Gram_Options opt;
  opt.prime_group = 32;
  opt.timing = flag(argc, argv, "--timing");
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
  printf("device %s; Block_Factor batch=%d n=%d, B: n x %d, column scale 2^-%.2f j, "
         "refine %d\n", prop.name, batch, n, nc, cs, refine);
  Block_Factor F(n, batch, L, D, opt);
  Mpf_Batch Mm = F.alloc(n, n), S, Bm = F.alloc(n, nc), Lf = F.alloc(n, n),
            X = F.alloc(n, n), P, Id = F.alloc(n, n);
  fill_rand<<<(Mm.entries() + 255) / 256, 256>>>(Mm.m, Mm.e, Mm.entries(), L, n,
                                                 n, 7, cs, row_span);
  fill_rand<<<(Bm.entries() + 255) / 256, 256>>>(Bm.m, Bm.e, Bm.entries(), L, n,
                                                 nc, 8, 0.0, 0.0);
  F.mm(Mm, false, Mm, true, S); // S = M M^T
  // + 2^-40 n I: keeps S SPD after rounding without changing its conditioning much
  fill_diag<<<(Id.entries() + 255) / 256, 256>>>(Id.m, Id.e, Id.entries(), L, n,
                                                 -ridge + int64_t(std::log2(n)));
  mpf_batch_add(S, Id, 1, S, 0);
  F.release(Mm);
  F.release(Id);
  F.clear_cache();
  Mpf_Batch Sgpu = S;
  if(perturb > 0)
    {
      Mpf_Batch E = F.alloc(n, n);
      fill_diag<<<(E.entries() + 255) / 256, 256>>>(E.m, E.e, E.entries(), L, n,
                                                   -int64_t(perturb));
      Sgpu = F.alloc(n, n);
      mpf_batch_add(S, E, 1, Sgpu, 0);
      F.release(E);
      printf("negative control: the GPU factors S + 2^-%d I; the checks use S\n",
             perturb);
    }
  cudaError_t ce = cudaDeviceSynchronize();
  if(ce != cudaSuccess)
    {
      printf("CUDA error %s\n", cudaGetErrorString(ce));
      return 2;
    }
  double bestF = 1e30, bestP = 1e30, bestR = 1e30;
  for(int rep = 0; rep < reps; ++rep)
    {
      cudaDeviceSynchronize();
      auto t0 = std::chrono::steady_clock::now();
      F.factor(Sgpu, Lf, X, fo);
      cudaDeviceSynchronize();
      auto t1 = std::chrono::steady_clock::now();
      if(drop)
        F.clear_cache();
      if(P.m)
        F.release(P);
      F.solve(Lf, X, Bm, P, 0);
      cudaDeviceSynchronize();
      auto t2 = std::chrono::steady_clock::now();
      if(refine)
        {
          F.release(P);
          F.solve(Lf, X, Bm, P, refine);
          cudaDeviceSynchronize();
        }
      auto t3 = std::chrono::steady_clock::now();
      double f = std::chrono::duration<double>(t1 - t0).count(),
             p = std::chrono::duration<double>(t2 - t1).count(),
             r = std::chrono::duration<double>(t3 - t2).count();
      printf("rep %d: factor %.4f s, P = X B %.4f s, P with %d refinement(s) %.4f s\n",
             rep, f, p, refine, r);
      bestF = std::min(bestF, f);
      bestP = std::min(bestP, p);
      bestR = std::min(bestR, r);
      ce = cudaDeviceSynchronize();
      if(ce != cudaSuccess)
        {
          printf("CUDA error %s\n", cudaGetErrorString(ce));
          return 2;
        }
    }
  printf("best: factor %.4f s, X B %.4f s, refined solve %.4f s (first rep includes "
         "allocation)\n", bestF, bestP, bestR);
  for(auto &[k, v] : F.timings())
    printf("  phase %-28s %9.2f ms (all reps)\n", k.c_str(), v);

  // ---- checks
  Mpf_Batch Zq;
  if(solve_q)
    {
      F.mm(X, true, P, false, Zq); // X^T P = L^-T L^-1 B
      F.check();
      cudaDeviceSynchronize();
    }
  mpf_set_default_prec(2 * 32 * L + 256);
  Host hS = down(S), hL = down(Lf), hX = down(X), hB = down(Bm), hP = down(P);
  Host hZ = solve_q ? down(Zq) : Host{};
  const int64_t nn = int64_t(n) * n, nbk = int64_t(n) * nc;
  std::mt19937_64 rng(5);
  mpf_t a, b, c, s, t, mx;
  mpf_inits(a, b, c, s, t, mx, nullptr);
  double res_S = -1e9, res_P = -1e9, res_X = -1e9;
  int nchk = checks > 0 ? checks : 200;
  for(int q = 0; q < nchk; ++q)
    {
      int64_t bb = rng() % batch, i = rng() % n, j = rng() % n;
      if(i < j)
        std::swap(i, j);
      // (L L^T)_ij - S_ij
      mpf_set_ui(s, 0);
      for(int64_t k = 0; k <= j; ++k)
        {
          val(a, hL, bb * nn + k * n + i);
          val(b, hL, bb * nn + k * n + j);
          mpf_mul(t, a, b);
          mpf_add(s, s, t);
        }
      val(a, hS, bb * nn + j * n + i);
      mpf_sub(s, s, a);
      val(c, hS, bb * nn + i * n + i);
      val(t, hS, bb * nn + j * n + j);
      mpf_mul(c, c, t);
      mpf_sqrt(c, c); // sqrt(S_ii S_jj)
      res_S = std::max(res_S, lg(s) - lg(c));
      // (X L - I)_ij relative to max_k |X_ik| |L_kj| -> use |X_ii L_ii| = 1 scale
      mpf_set_ui(s, 0);
      for(int64_t k = j; k <= i; ++k)
        {
          val(a, hX, bb * nn + k * n + i);
          val(b, hL, bb * nn + j * n + k);
          mpf_mul(t, a, b);
          mpf_add(s, s, t);
        }
      if(i == j)
        mpf_sub_ui(s, s, 1);
      res_X = std::max(res_X, lg(s));
      // (L P - B)_ic relative to |L_i:| |P_:c|
      int64_t col = rng() % nc;
      mpf_set_ui(s, 0);
      double lmax = -1e9, pmax = -1e9;
      for(int64_t k = 0; k <= i; ++k)
        {
          val(a, hL, bb * nn + k * n + i);
          val(b, hP, bb * nbk + col * n + k);
          lmax = std::max(lmax, lg(a));
          pmax = std::max(pmax, lg(b));
          mpf_mul(t, a, b);
          mpf_add(s, s, t);
        }
      val(a, hB, bb * nbk + col * n + i);
      mpf_sub(s, s, a);
      res_P = std::max(res_P, lg(s) - lmax - pmax - std::log2(double(i + 1)));
    }
  printf("residuals: max|LL^T - S|/sqrt(S_ii S_jj) = 2^%.1f, max|XL - I| = 2^%.1f, "
         "max|LP - B|/(k max|L| max|P|) = 2^%.1f\n", res_S, res_X, res_P);
  // -850 at L = 28
  bool ok = res_S < -(32.0 * L - 46) && res_P < -(32.0 * L - 46);
  if(ref)
    {
      // GMP Cholesky of block 0 and substitution for all columns
      std::vector<__mpf_struct> R(nn);
      for(int64_t q = 0; q < nn; ++q)
        {
          mpf_init(&R[q]);
          val(&R[q], hS, q);
        }
      for(int64_t j = 0; j < n; ++j)
        {
          mpf_set(s, &R[j * n + j]);
          for(int64_t k = 0; k < j; ++k)
            {
              mpf_mul(t, &R[k * n + j], &R[k * n + j]);
              mpf_sub(s, s, t);
            }
          mpf_sqrt(&R[j * n + j], s);
          for(int64_t i = j + 1; i < n; ++i)
            {
              mpf_set(s, &R[j * n + i]);
              for(int64_t k = 0; k < j; ++k)
                {
                  mpf_mul(t, &R[k * n + i], &R[k * n + j]);
                  mpf_sub(s, s, t);
                }
              mpf_div(&R[j * n + i], s, &R[j * n + j]);
            }
        }
      double fwL = -1e9, fwP = -1e9, fwZ = -1e9, lref = -1e9;
      for(int64_t j = 0; j < n; ++j)
        for(int64_t i = j; i < n; ++i)
          lref = std::max(lref, lg(&R[j * n + i]));
      for(int64_t j = 0; j < n; ++j)
        for(int64_t i = j; i < n; ++i)
          {
            val(a, hL, j * n + i);
            mpf_sub(a, a, &R[j * n + i]);
            fwL = std::max(fwL, lg(a) - lref);
          }
      std::vector<__mpf_struct> y(n);
      for(int64_t i = 0; i < n; ++i)
        mpf_init(&y[i]);
      for(int64_t col = 0; col < std::min(nc, 8); ++col)
        {
          double pref = -1e9;
          for(int64_t i = 0; i < n; ++i)
            {
              val(s, hB, col * n + i);
              for(int64_t k = 0; k < i; ++k)
                {
                  mpf_mul(t, &R[k * n + i], &y[k]);
                  mpf_sub(s, s, t);
                }
              mpf_div(&y[i], s, &R[i * n + i]);
              pref = std::max(pref, lg(&y[i]));
            }
          for(int64_t i = 0; i < n; ++i)
            {
              val(a, hP, col * n + i);
              mpf_sub(a, a, &y[i]);
              fwP = std::max(fwP, lg(a) - pref);
            }
          if(solve_q)
            {
              // z = L^-T y by backward substitution (R[k*n + i] = L_ik)
              double zref = -1e9;
              for(int64_t i = n - 1; i >= 0; --i)
                {
                  mpf_set(s, &y[i]);
                  for(int64_t k = i + 1; k < n; ++k)
                    {
                      mpf_mul(t, &R[i * n + k], &y[k]);
                      mpf_sub(s, s, t);
                    }
                  mpf_div(&y[i], s, &R[i * n + i]); // y now holds z for k >= i
                  zref = std::max(zref, lg(&y[i]));
                }
              for(int64_t i = 0; i < n; ++i)
                {
                  val(a, hZ, col * n + i);
                  mpf_sub(a, a, &y[i]);
                  fwZ = std::max(fwZ, lg(a) - zref);
                }
            }
        }
      ok = ok && fwL < forward_limit && fwP < forward_limit && (!solve_q || fwZ < forward_limit);
      if(solve_q)
        printf("solve with S through the inverse (block 0, vs GMP): max|X^T X B - S^-1 B|/max|S^-1 B| = 2^%.1f "
               "(per column, up to 8 columns)\n", fwZ);
      printf("forward (block 0, vs GMP at %d bits): max|L - L_ref|/max|L_ref| = 2^%.1f, "
             "max|P - P_ref|/max|P_ref| = 2^%.1f\n", 2 * 32 * L + 256, fwL, fwP);
    }
  printf("CHECK %s\n", ok ? "PASSED" : "FAILED");
  return ok ? 0 : 1;
}
