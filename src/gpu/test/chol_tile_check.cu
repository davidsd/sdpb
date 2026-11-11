// GPU backend for SDPB: check of the ported base case (Tile_Kernels: batched
// Cholesky and triangular inverse of 32 x 32 tiles in wide fixed point).
//
// Random SPD tiles A = M M^T with column scales 2^(-c j) (condition number up
// to ~2^(2c*31)), equilibrated to |a_ij| <= 1, rounded to F = 32 LF fractional
// bits.  Checks, all in exact GMP integer arithmetic on the host:
//   backward error   max |L L^T - A|                (expect ~ BT 2^-F)
//   inverse residual max |X L - I| / max|X|          (expect ~ BT 2^-F)
//   forward error    max |L - L_ref| against a GMP mpf Cholesky at 2F+128 bits
//                    of the same A (independent of the GPU path)
//   chol_tile_check --tiles 64 --scale 3 [--fused] [--reps 3] [--prof]
// --fused: the fused leaf (leaf_fused.cuh) instead of chol_tiles +
// trinv_tiles; the checks are the same.

#include "../rns/Limb_Sizes.hxx"
#include "../linalg/basecase.cuh"
#include "../linalg/leaf_fused.cuh"

#include <gmp.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

using namespace sdpb_gpu;

// make CHECK_LF=<L>: one tile size per binary
#ifndef CHECK_LF
#define CHECK_LF 28
#endif
constexpr int LF = CHECK_LF, WL = leaf_limbs(LF), BT = 32;

namespace
{
  void to_limbs(const mpz_t x, uint32_t *w)
  {
    mpz_t t, m;
    mpz_inits(t, m, nullptr);
    mpz_ui_pow_ui(m, 2, 32 * WL);
    mpz_mod(t, x, m); // two's complement
    for(int i = 0; i < WL; ++i)
      w[i] = 0;
    size_t cnt = 0;
    mpz_export(w, &cnt, -1, 4, 0, 0, t);
    mpz_clears(t, m, nullptr);
  }
  void from_limbs(mpz_t x, const uint32_t *w)
  {
    mpz_import(x, WL, -1, 4, 0, 0, w);
    if(w[WL - 1] >> 31)
      {
        mpz_t m;
        mpz_init(m);
        mpz_ui_pow_ui(m, 2, 32 * WL);
        mpz_sub(x, x, m);
        mpz_clear(m);
      }
  }
  double lg(const mpz_t x)
  {
    if(mpz_sgn(x) == 0)
      return -1e9;
    long e;
    double d = mpz_get_d_2exp(&e, x);
    return std::log2(std::fabs(d)) + e;
  }
  const char *arg(int argc, char **argv, const char *name, const char *dflt)
  {
    for(int i = 1; i + 1 < argc; ++i)
      if(!strcmp(argv[i], name))
        return argv[i + 1];
    return dflt;
  }
}

int main(int argc, char **argv)
{
  const int nt = atoi(arg(argc, argv, "--tiles", "64"));
  const double c = atof(arg(argc, argv, "--scale", "3"));
  bool fused = false, prof = false, thread_leaf = false, no_la = false;
  for(int i = 1; i < argc; ++i)
    {
      fused = fused || !strcmp(argv[i], "--fused");
      prof = prof || !strcmp(argv[i], "--prof");
      thread_leaf = thread_leaf || !strcmp(argv[i], "--thread-leaf");
      no_la = no_la || !strcmp(argv[i], "--no-lookahead");
    }
  const int reps = atoi(arg(argc, argv, "--reps", "1"));
  constexpr int NT = 512;
  const int F = 32 * LF;
  std::mt19937_64 rng(42);
  std::vector<uint32_t> hA(size_t(nt) * BT * BT * WL);
  std::vector<__mpz_struct> A(size_t(nt) * BT * BT);
  mpz_t t, u, v;
  mpz_inits(t, u, v, nullptr);
  for(int q = 0; q < nt; ++q)
    {
      // M with 2^60-scale random integers, columns scaled by 2^(-c j)
      std::vector<__mpz_struct> M(BT * BT);
      for(int i = 0; i < BT; ++i)
        for(int j = 0; j < BT; ++j)
          {
            mpz_init(&M[i * BT + j]);
            int64_t r = int64_t(rng() >> 4) - (int64_t(1) << 59);
            mpz_set_si(&M[i * BT + j], r);
            mpz_mul_2exp(&M[i * BT + j], &M[i * BT + j],
                         (unsigned long)std::lround(c * (BT - 1 - j)));
          }
      // A = M M^T, then scale so max |a_ij| <= 1 at F fractional bits
      std::vector<__mpz_struct> G(BT * BT);
      double top = -1e9;
      for(int i = 0; i < BT; ++i)
        for(int j = 0; j < BT; ++j)
          {
            mpz_init(&G[i * BT + j]);
            for(int k = 0; k < BT; ++k)
              mpz_addmul(&G[i * BT + j], &M[i * BT + k], &M[j * BT + k]);
            top = std::max(top, lg(&G[i * BT + j]));
          }
      long sh = long(std::ceil(top)) + 1 - F; // a = G 2^-sh, |a| < 2^F
      for(int e = 0; e < BT * BT; ++e)
        {
          __mpz_struct *a = &A[size_t(q) * BT * BT + e];
          mpz_init(a);
          if(sh > 0)
            mpz_fdiv_q_2exp(a, &G[e], sh);
          else
            mpz_mul_2exp(a, &G[e], -sh);
          to_limbs(a, &hA[(size_t(q) * BT * BT + e) * WL]);
          mpz_clear(&G[e]);
          mpz_clear(&M[e]);
        }
    }
  uint32_t *dA, *dA0, *dR, *dX, *dS;
  long long *dP;
  size_t bytes = hA.size() * 4;
  cudaMalloc(&dA, bytes);
  cudaMalloc(&dA0, bytes);
  cudaMalloc(&dX, bytes);
  cudaMalloc(&dS, bytes);
  cudaMalloc(&dR, size_t(nt) * BT * WL * 4);
  cudaMalloc(&dP, size_t(nt) * 3 * 8);
  cudaMemcpy(dA0, hA.data(), bytes, cudaMemcpyHostToDevice);
  cudaEvent_t e0, e1, e2;
  cudaEventCreate(&e0);
  cudaEventCreate(&e1);
  cudaEventCreate(&e2);
  float ms1 = 1e30f, ms2 = 1e30f;
  for(int r = 0; r < reps; ++r)
    {
      cudaMemcpy(dA, dA0, bytes, cudaMemcpyDeviceToDevice);
      cudaEventRecord(e0);
      if(fused)
        {
          if(thread_leaf) // pivot and row/column k one thread each
            leaf_chol_inv<WL, LF, BT, NT, false><<<nt, NT>>>(dA, dX, prof ? dP : nullptr);
          else if(no_la) // pivot between the phases
            leaf_chol_inv<WL, LF, BT, NT, true, false><<<nt, NT>>>(dA, dX, prof ? dP : nullptr);
          else
            leaf_chol_inv<WL, LF, BT, NT, true, true><<<nt, NT>>>(dA, dX, prof ? dP : nullptr);
          cudaEventRecord(e1);
        }
      else
        {
          chol_tiles<WL, LF, BT><<<nt, 256>>>(dA, dR);
          cudaEventRecord(e1);
          // trinv reads only i >= j
          trinv_tiles<WL, LF, BT><<<nt, 256>>>(dA, dR, dX, dS);
        }
      cudaEventRecord(e2);
      cudaError_t err = cudaDeviceSynchronize();
      if(err != cudaSuccess)
        {
          printf("CUDA error %s\n", cudaGetErrorString(err));
          return 2;
        }
      float a, b;
      cudaEventElapsedTime(&a, e0, e1);
      cudaEventElapsedTime(&b, e1, e2);
      ms1 = std::min(ms1, a);
      ms2 = std::min(ms2, b);
    }
  printf("tiles %d of %dx%d, WL=%d LF=%d (F=%d bits), column scale 2^-%.1f j, "
         "%s: chol %.3f ms, inverse %.3f ms, total %.3f ms (best of %d)\n",
         nt, BT, BT, WL, LF, F, c, fused ? "fused leaf" : "chol_tiles + trinv_tiles",
         ms1, ms2, ms1 + ms2, reps);
  if(fused && prof)
    {
      std::vector<long long> hp(size_t(nt) * 3);
      cudaMemcpy(hp.data(), dP, hp.size() * 8, cudaMemcpyDeviceToHost);
      printf("tile 0 cycles: A (pivot) %lld, B (row/column) %lld, C (updates) %lld\n",
             hp[0], hp[1], hp[2]);
    }
  std::vector<uint32_t> hL(hA.size()), hX(hA.size());
  cudaMemcpy(hL.data(), dA, bytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(hX.data(), dX, bytes, cudaMemcpyDeviceToHost);

  double worst_bw = -1e9, worst_inv = -1e9, worst_fw = -1e9, worst_cond = -1e9;
  mpf_set_default_prec(2 * F + 128);
  for(int q = 0; q < nt; ++q)
    {
      std::vector<__mpz_struct> Lz(BT * BT), Xz(BT * BT);
      for(int e = 0; e < BT * BT; ++e)
        {
          mpz_init(&Lz[e]);
          mpz_init(&Xz[e]);
          int i = e / BT, j = e % BT;
          if(i >= j)
            {
              from_limbs(&Lz[e], &hL[(size_t(q) * BT * BT + e) * WL]);
              from_limbs(&Xz[e], &hX[(size_t(q) * BT * BT + e) * WL]);
            }
        }
      // backward error: (L L^T)_ij 2^-2F - a_ij 2^-F
      for(int i = 0; i < BT; ++i)
        for(int j = 0; j <= i; ++j)
          {
            mpz_set_ui(t, 0);
            for(int k = 0; k <= j; ++k)
              mpz_addmul(t, &Lz[i * BT + k], &Lz[j * BT + k]);
            mpz_mul_2exp(u, &A[size_t(q) * BT * BT + i * BT + j], F);
            mpz_sub(t, t, u);
            worst_bw = std::max(worst_bw, lg(t) - 2.0 * F);
          }
      // inverse residual (X L - I) relative to max |X|
      double xmax = -1e9;
      for(int e = 0; e < BT * BT; ++e)
        xmax = std::max(xmax, lg(&Xz[e]) - F);
      worst_cond = std::max(worst_cond, xmax);
      for(int i = 0; i < BT; ++i)
        for(int j = 0; j <= i; ++j)
          {
            mpz_set_ui(t, 0);
            for(int k = j; k <= i; ++k)
              mpz_addmul(t, &Xz[i * BT + k], &Lz[k * BT + j]);
            if(i == j)
              {
                mpz_ui_pow_ui(u, 2, 2 * F);
                mpz_sub(t, t, u);
              }
            worst_inv = std::max(worst_inv, lg(t) - 2.0 * F - xmax);
          }
      // forward error against a GMP mpf Cholesky of the same A
      std::vector<__mpf_struct> Rv(BT * BT);
      for(int e = 0; e < BT * BT; ++e)
        {
          mpf_init((&Rv[e]));
          mpf_set_z((&Rv[e]), &A[size_t(q) * BT * BT + e]);
          mpf_div_2exp((&Rv[e]), (&Rv[e]), F);
        }
      mpf_t s, w;
      mpf_inits(s, w, nullptr);
      for(int j = 0; j < BT; ++j)
        {
          mpf_set(s, (&Rv[j * BT + j]));
          for(int k = 0; k < j; ++k)
            {
              mpf_mul(w, (&Rv[j * BT + k]), (&Rv[j * BT + k]));
              mpf_sub(s, s, w);
            }
          mpf_sqrt((&Rv[j * BT + j]), s);
          for(int i = j + 1; i < BT; ++i)
            {
              mpf_set(s, (&Rv[i * BT + j]));
              for(int k = 0; k < j; ++k)
                {
                  mpf_mul(w, (&Rv[i * BT + k]), (&Rv[j * BT + k]));
                  mpf_sub(s, s, w);
                }
              mpf_div((&Rv[i * BT + j]), s, (&Rv[j * BT + j]));
            }
        }
      for(int i = 0; i < BT; ++i)
        for(int j = 0; j <= i; ++j)
          {
            mpf_set_z(w, &Lz[i * BT + j]);
            mpf_div_2exp(w, w, F);
            mpf_sub(w, w, (&Rv[i * BT + j]));
            mpf_abs(w, w);
            if(mpf_sgn(w) != 0)
              {
                long ex;
                double d = mpf_get_d_2exp(&ex, w);
                worst_fw = std::max(worst_fw, std::log2(d) + ex);
              }
          }
      for(int e = 0; e < BT * BT; ++e)
        {
          mpf_clear((&Rv[e]));
          mpz_clear(&Lz[e]);
          mpz_clear(&Xz[e]);
        }
      mpf_clears(s, w, nullptr);
    }
  bool ok = worst_bw < -(F - 8) && worst_inv < -(F - 8 - worst_cond * 0)
            && worst_fw < -(F - 8) + 2 * worst_cond;
  printf("CHECK %s: max|LL^T - A| = 2^%.1f, max|XL - I|/max|X| = 2^%.1f, "
         "max|L - L_ref| = 2^%.1f (F = %d; max|X| = 2^%.1f)\n",
         ok ? "PASSED" : "FAILED", worst_bw, worst_inv, worst_fw, F,
         worst_cond);
  return ok ? 0 : 1;
}
