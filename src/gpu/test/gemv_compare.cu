// mpf_gemv on real blocks: one thread per
// output against several (Zechuan Zheng's partial sums).  Same job, same data:
// the shapes of an iteration, timed; each variant must be bitwise equal to one
// thread per output, and sampled outputs must agree with an exact GMP dot product
// to the bound relative to the largest term.
//   gemv_compare --dir <Schur dump> --nb 36 --P 945 --N 3406 --q-dir <Q dump>
//     [--targets 16384,65536,262144] [--reps 6] [--perturb-bits b]
// --perturb-bits b: GMP sees 2^-b added to the matrix entry of a term (negative
// control: the exact check must fail for b well above the bound).
#include "../linalg/Block_Factor.hxx"
#include "../solver/Mpf_Ops.hxx"

#include <gmp.h>

#include <algorithm>
#include <chrono>
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
  constexpr int L = 28, NB = 32 * L - 2;
  const char *arg(int argc, char **argv, const char *name, const char *dflt)
  {
    for(int i = 1; i + 1 < argc; ++i)
      if(!strcmp(argv[i], name))
        return argv[i + 1];
    return dflt;
  }
  template <class T> std::vector<T> slurp(const std::string &path, size_t count)
  {
    std::vector<T> v(count);
    FILE *f = fopen(path.c_str(), "rb");
    if(!f || fread(v.data(), sizeof(T), count, f) != count)
      throw std::runtime_error("read " + path);
    fclose(f);
    return v;
  }
  void check(cudaError_t e)
  {
    if(e != cudaSuccess)
      throw std::runtime_error(cudaGetErrorString(e));
  }
  struct Mat
  {
    Mpf_Batch b;
    Mat(int batch, int64_t r, int64_t c)
    {
      b.L = L;
      b.batch = batch;
      b.rows = r;
      b.cols = c;
      check(cudaMalloc(&b.m, size_t(L) * b.entries() * 4));
      check(cudaMalloc(&b.e, size_t(b.entries()) * 8));
    }
    ~Mat()
    {
      cudaFree(b.m);
      cudaFree(b.e);
    }
    Mat(const Mat &) = delete;
  };
  void up(Mat &M, const int32_t *m, const int64_t *e)
  {
    check(cudaMemcpy(M.b.m, m, size_t(L) * M.b.entries() * 4, cudaMemcpyHostToDevice));
    check(cudaMemcpy(M.b.e, e, size_t(M.b.entries()) * 8, cudaMemcpyHostToDevice));
  }
  struct Host
  {
    std::vector<int32_t> m;
    std::vector<int64_t> e;
    int64_t n;
    void get(const Mpf_Batch &A)
    {
      n = A.entries();
      m.resize(size_t(L) * n);
      e.resize(n);
      check(cudaMemcpy(m.data(), A.m, size_t(L) * n * 4, cudaMemcpyDeviceToHost));
      check(cudaMemcpy(e.data(), A.e, size_t(n) * 8, cudaMemcpyDeviceToHost));
    }
  };
  void to_mpz(mpz_t z, const Host &H, int64_t i)
  {
    uint32_t w[L];
    for(int t = 0; t < L; ++t)
      w[t] = uint32_t(H.m[int64_t(t) * H.n + i]);
    mpz_import(z, L, -1, 4, 0, 0, w);
    if(w[L - 1] >> 31)
      {
        mpz_t p;
        mpz_init(p);
        mpz_setbit(p, 32 * L);
        mpz_sub(z, z, p);
        mpz_clear(p);
      }
  }
  double gemv_log2_err(const Host &A, const Host &v, const Host &O, const Gemv_Spec &g,
                       int64_t b, int64_t o, int perturb_bits)
  {
    // exact sum of the terms as an integer times 2^emin; error relative to the largest term
    const int64_t a0 = b * g.ab + o * g.ao, v0 = b * g.vb + o * g.vo;
    int64_t top = INT64_MIN;
    std::vector<std::pair<int64_t, int64_t>> terms;
    for(int64_t tb = 0; tb < g.ntb; ++tb)
      for(int64_t ti = 0; ti < g.nti; ++ti)
        {
          const int64_t ia = a0 + tb * g.atb + ti * g.ati, iv = v0 + tb * g.vtb + ti * g.vti;
          if(A.e[ia] == MPF_EZ || v.e[iv] == MPF_EZ)
            continue;
          terms.emplace_back(ia, iv);
          top = std::max(top, A.e[ia] + v.e[iv]);
        }
    mpz_t sum, x, y, t, r;
    mpz_inits(sum, x, y, t, r, nullptr);
    const int64_t base = top - 3 * NB; // terms far below contribute under the bound anyway
    for(auto &tm : terms)
      {
        to_mpz(x, A, tm.first);
        to_mpz(y, v, tm.second);
        if(perturb_bits > 0 && A.e[tm.first] + v.e[tm.second] == top)
          {
            mpz_set(t, x);
            mpz_tdiv_q_2exp(t, t, perturb_bits);
            mpz_add(x, x, t);
          }
        mpz_mul(t, x, y);
        const int64_t sh = A.e[tm.first] + v.e[tm.second] - base;
        if(sh >= 0)
          mpz_mul_2exp(t, t, sh);
        else
          mpz_tdiv_q_2exp(t, t, -sh);
        mpz_add(sum, sum, t);
      }
    const int64_t iout = b * g.ob + o * g.oo;
    if(O.e[iout] == MPF_EZ)
      mpz_set_ui(r, 0);
    else
      {
        to_mpz(r, O, iout);
        const int64_t sh = O.e[iout] - base;
        if(sh >= 0)
          mpz_mul_2exp(r, r, sh);
        else
          mpz_tdiv_q_2exp(r, r, -sh);
      }
    mpz_sub(t, sum, r);
    long ex = 0;
    const double d = mpz_sgn(t) == 0 ? 0 : mpz_get_d_2exp(&ex, t);
    mpz_clears(sum, x, y, t, r, nullptr);
    if(d == 0)
      return -100000;
    // log2 |error| at 2^base, against the largest term's scale 2^(top + 2 NB)
    return ex + std::log2(std::fabs(d)) + double(base) - double(top + 2 * NB);
  }
}

int main(int argc, char **argv)
{
  const std::string dir = arg(argc, argv, "--dir", "."), qdir = arg(argc, argv, "--q-dir", "");
  const int nb = atoi(arg(argc, argv, "--nb", "36")), P = atoi(arg(argc, argv, "--P", "945")),
            N = atoi(arg(argc, argv, "--N", "3406")), reps = atoi(arg(argc, argv, "--reps", "6"));
  const int Qn = atoi(arg(argc, argv, "--q-n", "3406"));
  const int perturb = atoi(arg(argc, argv, "--perturb-bits", "0"));
  std::vector<int> targets;
  {
    std::string s = arg(argc, argv, "--targets", "65536,262144");
    for(size_t i = 0; i < s.size();)
      {
        size_t j = s.find(',', i);
        if(j == std::string::npos)
          j = s.size();
        targets.push_back(atoi(s.substr(i, j - i).c_str()));
        i = j + 1;
      }
  }
  const int64_t R = int64_t(nb) * P, bp = R * N;
  Mat B(nb, P, N), Pf(1, R, N);
  {
    // the dumps are row-major (SDPviBe); the backend is column-major
    auto m = slurp<int32_t>(dir + "/B0.m.bin", size_t(L) * bp);
    auto e = slurp<int64_t>(dir + "/B0.e.bin", bp);
    std::vector<int32_t> m2(size_t(L) * bp);
    std::vector<int64_t> e2(bp);
    for(int b = 0; b < nb; ++b)
      for(int i = 0; i < P; ++i)
        for(int j = 0; j < N; ++j)
          {
            const int64_t src = (int64_t(b) * P + i) * N + j, dst = int64_t(b) * P * N + int64_t(j) * P + i;
            e2[dst] = e[src];
            for(int t = 0; t < L; ++t)
              m2[int64_t(t) * bp + dst] = m[int64_t(t) * bp + src];
          }
    up(B, m2.data(), e2.data());
    m = slurp<int32_t>(dir + "/P0.m.bin", size_t(L) * bp);
    e = slurp<int64_t>(dir + "/P0.e.bin", bp);
    for(int64_t r = 0; r < R; ++r)
      for(int j = 0; j < N; ++j)
        {
          const int64_t src = r * N + j, dst = int64_t(j) * R + r;
          e2[dst] = e[src];
          for(int t = 0; t < L; ++t)
            m2[int64_t(t) * bp + dst] = m[int64_t(t) * bp + src];
        }
    up(Pf, m2.data(), e2.data());
  }
  // real factors: the Schur blocks' and Q's inverses, from Block_Factor
  const int np = (P + 31) / 32 * 32, qp = (Qn + 31) / 32 * 32;
  const int64_t sp = int64_t(nb) * np * np;
  Mat Xs(nb, np, np), Xq(1, qp, qp);
  {
    auto m = slurp<int32_t>(dir + "/S0.m.bin", size_t(L) * nb * int64_t(P) * P);
    auto e = slurp<int64_t>(dir + "/S0.e.bin", int64_t(nb) * P * P);
    Mat S(nb, P, P), Sp(nb, np, np), Lf(nb, np, np);
    std::vector<int32_t> m2(size_t(L) * nb * int64_t(P) * P);
    std::vector<int64_t> e2(int64_t(nb) * P * P);
    const int64_t tot = int64_t(nb) * P * P;
    for(int b = 0; b < nb; ++b)
      for(int i = 0; i < P; ++i)
        for(int j = 0; j < P; ++j)
          {
            const int64_t src = int64_t(b) * P * P + int64_t(i) * P + j, dst = int64_t(b) * P * P + int64_t(j) * P + i;
            e2[dst] = e[src];
            for(int t = 0; t < L; ++t)
              m2[int64_t(t) * tot + dst] = m[int64_t(t) * tot + src];
          }
    up(S, m2.data(), e2.data());
    mpf_pad_identity(S.b, Sp.b, 0);
    Block_Factor F(32, nb, L, 112);
    F.factor(Sp.b, Lf.b, Xs.b);
    F.check();
    check(cudaStreamSynchronize(0));
  }
  {
    auto m = slurp<int32_t>(qdir + "/S0.m.bin", size_t(L) * int64_t(Qn) * Qn);
    auto e = slurp<int64_t>(qdir + "/S0.e.bin", int64_t(Qn) * Qn);
    Mat Q(1, Qn, Qn), Qp(1, qp, qp), Lf(1, qp, qp);
    std::vector<int32_t> m2(size_t(L) * int64_t(Qn) * Qn);
    std::vector<int64_t> e2(int64_t(Qn) * Qn);
    const int64_t tot = int64_t(Qn) * Qn;
    for(int i = 0; i < Qn; ++i)
      for(int j = 0; j < Qn; ++j)
        {
          const int64_t src = int64_t(i) * Qn + j, dst = int64_t(j) * Qn + i;
          e2[dst] = e[src];
          for(int t = 0; t < L; ++t)
            m2[int64_t(t) * tot + dst] = m[int64_t(t) * tot + src];
        }
    up(Q, m2.data(), e2.data());
    mpf_pad_identity(Q.b, Qp.b, 0);
    Block_Factor F(32, 1, L, 112);
    F.factor(Qp.b, Lf.b, Xq.b, Factor_Options{true, 1, 0, -1, 1024, true});
    F.check();
    check(cudaStreamSynchronize(0));
  }
  // Vectors of real data, dense as w and dy are: entry r of the R-vector is the
  // first nonzero entry of row r of P from column N / 3 on, entry c of the N-vector
  // the first nonzero entry of column c from row P / 2 on.  (One column or row of P
  // alone is mostly exact zeros, which would time little but exponent reads.)
  Host HP, HB, HXs, HXq;
  HP.get(Pf.b);
  HB.get(B.b);
  HXs.get(Xs.b);
  HXq.get(Xq.b);
  Mat wR(1, R, 1), wB(nb, P, 1), vN(1, N, 1), wS(nb, np, 1), vQ(1, qp, 1);
  {
    std::vector<int32_t> m(size_t(L) * R, 0);
    std::vector<int64_t> e(R, MPF_EZ);
    int64_t filled = 0;
    for(int64_t r = 0; r < R; ++r)
      for(int d = 0; d < N; ++d)
        {
          const int64_t q = int64_t((N / 3 + d) % N) * R + r;
          if(HP.e[q] == MPF_EZ)
            continue;
          e[r] = HP.e[q];
          for(int t = 0; t < L; ++t)
            m[int64_t(t) * R + r] = HP.m[int64_t(t) * HP.n + q];
          ++filled;
          break;
        }
    printf("VECTORS w: %ld of %ld entries nonzero\n", long(filled), long(R));
    up(wR, m.data(), e.data());
    up(wB, m.data(), e.data());
    std::vector<int32_t> m3(size_t(L) * N);
    std::vector<int64_t> e3(N);
    filled = 0;
    for(int c = 0; c < N; ++c)
      {
        e3[c] = MPF_EZ;
        for(int64_t d = 0; d < R; ++d)
          {
            const int64_t q = int64_t(c) * R + (P / 2 + d) % R;
            if(HP.e[q] == MPF_EZ)
              continue;
            e3[c] = HP.e[q];
            for(int t = 0; t < L; ++t)
              m3[int64_t(t) * N + c] = HP.m[int64_t(t) * HP.n + q];
            ++filled;
            break;
          }
      }
    printf("VECTORS v: %ld of %d entries nonzero\n", long(filled), N);
    up(vN, m3.data(), e3.data());
    // padded rows of the square shapes get zeros
    const int64_t ws = int64_t(nb) * np;
    std::vector<int32_t> m4(size_t(L) * ws, 0);
    std::vector<int64_t> e4(ws, MPF_EZ);
    for(int b = 0; b < nb; ++b)
      for(int i = 0; i < P; ++i)
        {
          const int64_t dst = int64_t(b) * np + i, src = int64_t(b) * P + i;
          e4[dst] = e[src];
          for(int t = 0; t < L; ++t)
            m4[int64_t(t) * ws + dst] = m[int64_t(t) * R + src];
        }
    up(wS, m4.data(), e4.data());
    std::vector<int32_t> m5(size_t(L) * qp, 0);
    std::vector<int64_t> e5(qp, MPF_EZ);
    for(int c = 0; c < std::min(N, Qn); ++c)
      {
        e5[c] = e3[c];
        for(int t = 0; t < L; ++t)
          m5[int64_t(t) * qp + c] = m3[int64_t(t) * N + c];
      }
    up(vQ, m5.data(), e5.data());
  }
  struct Shape
  {
    const char *name;
    const Mpf_Batch *A, *v;
    const Host *HA;
    Gemv_Spec g;
    int outb;
    int64_t outr;
    int per_iter;
  };
  std::vector<Shape> shapes;
  auto sq = [&](int64_t n, int batch, bool tr) {
    Gemv_Spec g;
    g.batch = batch;
    g.no = n;
    g.ab = n * n;
    g.ao = tr ? n : 1;
    g.ati = tr ? 1 : n;
    g.nti = n;
    g.vb = n;
    g.vti = 1;
    g.ob = n;
    g.oo = 1;
    return g;
  };
  {
    Gemv_Spec g; // B y
    g.batch = nb; g.no = P; g.ab = int64_t(P) * N; g.ao = 1; g.ati = P; g.nti = N; g.vti = 1; g.ob = P; g.oo = 1;
    shapes.push_back({"B y", &B.b, &vN.b, &HB, g, nb, P, 1});
    Gemv_Spec h; // B^T x
    h.no = N; h.ao = P; h.atb = int64_t(P) * N; h.ati = 1; h.ntb = nb; h.nti = P; h.vtb = P; h.vti = 1; h.oo = 1;
    shapes.push_back({"B^T x", &B.b, &wB.b, &HB, h, 1, N, 1});
    Gemv_Spec i; // P^T w
    i.no = N; i.ao = R; i.ati = 1; i.nti = R; i.vti = 1; i.oo = 1;
    shapes.push_back({"P^T w", &Pf.b, &wR.b, &HP, i, 1, N, 2});
    Gemv_Spec j; // P dy
    j.no = R; j.ao = 1; j.ati = R; j.nti = N; j.vti = 1; j.oo = 1;
    shapes.push_back({"P dy", &Pf.b, &vN.b, &HP, j, nb, P, 2});
    shapes.push_back({"S^-1 w", &Xs.b, &wS.b, &HXs, sq(np, nb, false), nb, np, 2});
    shapes.push_back({"S^-T w", &Xs.b, &wS.b, &HXs, sq(np, nb, true), nb, np, 2});
    shapes.push_back({"Xq v", &Xq.b, &vQ.b, &HXq, sq(qp, 1, false), 1, qp, 2});
    shapes.push_back({"Xq^T v", &Xq.b, &vQ.b, &HXq, sq(qp, 1, true), 1, qp, 2});
  }
  std::mt19937_64 rng(336);
  double total_old = 0, total_new = 0;
  int failures = 0;
  for(auto &sh : shapes)
    {
      const int64_t nout = int64_t(sh.g.batch) * sh.g.no;
      Mat O0(sh.outb, sh.outr, 1), O1(sh.outb, sh.outr, 1);
      Host Hv;
      Hv.get(*sh.v);
      auto timed = [&](Mat &O, int target) {
        gemv_thread_target = target;
        std::vector<double> ts;
        for(int r = 0; r < reps; ++r)
          {
            const auto t0 = std::chrono::steady_clock::now();
            mpf_gemv(*sh.A, *sh.v, O.b, sh.g, 0);
            check(cudaStreamSynchronize(0));
            ts.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
          }
        ts.erase(ts.begin());
        std::sort(ts.begin(), ts.end());
        return ts[ts.size() / 2];
      };
      const double t0 = timed(O0, 0);
      Host H0;
      H0.get(O0.b);
      {
        // what the shape really does: a shape whose outputs have no terms times nothing
        int64_t nzo = 0, nza = 0, nzv = 0;
        for(int64_t i = 0; i < H0.n; ++i)
          nzo += H0.e[i] != MPF_EZ;
        for(int64_t i = 0; i < sh.HA->n; ++i)
          nza += sh.HA->e[i] != MPF_EZ;
        for(int64_t i = 0; i < Hv.n; ++i)
          nzv += Hv.e[i] != MPF_EZ;
        printf("NONZERO %-8s outputs %ld of %ld, matrix %ld of %ld, vector %ld of %ld\n", sh.name, long(nzo),
               long(H0.n), long(nza), long(sh.HA->n), long(nzv), long(Hv.n));
      }
      double best = t0;
      printf("GEMV %-8s outputs %ld terms %ld x%d/iter: one thread per output %.4f s", sh.name, long(nout),
             long(sh.g.ntb * sh.g.nti), sh.per_iter, t0);
      for(int target : targets)
        {
          const double t = timed(O1, target);
          Host H1;
          H1.get(O1.b);
          const bool same = H1.m == H0.m && H1.e == H0.e;
          if(!same)
            ++failures;
          printf(" | target %d: %.4f s (%.2fx%s)", target, t, t0 / t, same ? "" : " NOT BITWISE");
          if(same)
            best = std::min(best, t);
        }
      printf("\n");
      total_old += sh.per_iter * t0;
      total_new += sh.per_iter * best;
      // the exact check on outputs of the thread-per-output result, and of the last variant
      double worst = -1e9;
      Host H1;
      H1.get(O1.b);
      for(int k = 0; k < 24; ++k)
        {
          const int64_t q = rng() % nout, b = q / sh.g.no, o = q % sh.g.no;
          worst = std::max({worst, gemv_log2_err(*sh.HA, Hv, H1, sh.g, b, o, perturb)});
        }
      // each aligned term is rounded once, half an ulp below the largest term's
      // scale: at most nt 2^-(NB + 1) in all
      const double bound = -NB + std::ceil(std::log2(double(sh.g.ntb * sh.g.nti)));
      printf("GMPCHECK %-8s worst log2 error relative to the largest term %.1f over 24 outputs, bound %.0f: %s\n",
             sh.name, worst, bound, worst <= bound ? "pass" : "FAIL");
      if(worst > bound)
        ++failures;
      fflush(stdout);
    }
  printf("TOTAL over one chunk's iteration (shapes x calls): one thread per output %.3f s, best variant %.3f s (%.2fx); failures %d\n",
         total_old, total_new, total_old / total_new, failures);
  return failures ? 1 : 0;
}
