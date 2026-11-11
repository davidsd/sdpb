// All-entry GPU residual; independent GMP recomputation of
// worst positions and seeded samples. Reads exact stored Q and candidate L.
// Layout credited to SDPviBe (Antoine Tilloy, MIT); device product is backend.
#include "../linalg/Block_Factor.hxx"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <gmp.h>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
using namespace sdpb_gpu;
static const char *arg(int ac, char **av, const char *k, const char *d) {
  for (int i = 1; i + 1 < ac; i++)
    if (!strcmp(av[i], k))
      return av[i + 1];
  return d;
}
template <class T> std::vector<T> read(const std::string &p, size_t n) {
  std::vector<T> v(n);
  FILE *f = fopen(p.c_str(), "rb");
  if (!f || fread(v.data(), sizeof(T), n, f) != n || fgetc(f) != EOF)
    throw std::runtime_error("size/read " + p);
  fclose(f);
  return v;
}
struct H {
  int n, L;
  std::vector<int32_t> m;
  std::vector<int64_t> e;
};
static H readh(const std::string &p, int n, int L) {
  return {n, L, read<int32_t>(p + ".m.bin", size_t(L) * n * n),
          read<int64_t>(p + ".e.bin", size_t(n) * n)};
}
static void get(mpf_t a, const H &h, int i, int j, mpz_t z) {
  size_t idx = size_t(j) * h.n + i;
  if (h.e[idx] == MPF_EZ) {
    mpf_set_ui(a, 0);
    return;
  }
  uint32_t w[64];
  for (int t = 0; t < h.L; t++)
    w[t] = h.m[size_t(t) * h.n * h.n + idx];
  mpz_import(z, h.L, -1, 4, 0, 0, w);
  if (w[h.L - 1] >> 31) {
    mpz_t b;
    mpz_init_set_ui(b, 0);
    mpz_setbit(b, 32 * h.L);
    mpz_sub(z, z, b);
    mpz_clear(b);
  }
  mpf_set_z(a, z);
  if (h.e[idx] > 0)
    mpf_mul_2exp(a, a, h.e[idx]);
  else
    mpf_div_2exp(a, a, -h.e[idx]);
}
static double lg(mpf_t a) {
  if (!mpf_sgn(a))
    return -1e9;
  long e;
  double d = mpf_get_d_2exp(&e, a);
  return log2(fabs(d)) + e;
}
// Only a diagnostic: double-precision norms of stored multiprecision matrices.
// The all-entry product and independent GMP dots below check the inverse
// itself.
struct Norms {
  double one, inf, max_col_square, max_diag;
};
static Norms norms(const H &h, int P) {
  std::vector<double> rows(P, 0);
  Norms v{0, 0, 0, 0};
  const size_t plane = size_t(h.n) * h.n;
  for (int j = 0; j < P; j++) {
    double col = 0, sq = 0;
    for (int i = 0; i < P; i++) {
      size_t q = size_t(j) * h.n + i;
      if (h.e[q] == MPF_EZ)
        continue;
      double a = double(h.m[size_t(h.L - 1) * plane + q]) * 4294967296.0 +
                 double(uint32_t(h.m[size_t(h.L - 2) * plane + q]));
      a = std::abs(std::ldexp(a, int(h.e[q] + 32LL * (h.L - 2))));
      col += a;
      rows[i] += a;
      sq += a * a;
      if (i == j)
        v.max_diag = std::max(v.max_diag, a);
    }
    v.one = std::max(v.one, col);
    v.max_col_square = std::max(v.max_col_square, sq);
  }
  v.inf = *std::max_element(rows.begin(), rows.end());
  return v;
}
__device__ double logval(const int32_t *m, const int64_t *e, int L,
                         int64_t plane, int64_t q) {
  if (e[q] == MPF_EZ)
    return -1e9;
  double d = double(m[(L - 1) * plane + q]) * 4294967296.0 +
             double(uint32_t(m[(L - 2) * plane + q]));
  return log2(fabs(d)) + double(e[q] + 32LL * (L - 2));
}
__global__ void reduce_res(const int32_t *rm, const int64_t *re,
                           const int32_t *sm, const int64_t *se, int L,
                           int64_t n, int P, double *best, int64_t *pos) {
  __shared__ double v[256];
  __shared__ int64_t ix[256];
  int t = threadIdx.x;
  int64_t q = int64_t(blockIdx.x) * 256 + t, i = q % n, j = q / n;
  double x = -INFINITY;
  if (q < n * n && i < P && j <= i)
    x = logval(rm, re, L, n * n, q) -
        .5 * (logval(sm, se, L, n * n, i * n + i) +
              logval(sm, se, L, n * n, j * n + j));
  v[t] = x;
  ix[t] = q;
  __syncthreads();
  for (int d = 128; d; d /= 2) {
    if (t < d && v[t + d] > v[t]) {
      v[t] = v[t + d];
      ix[t] = ix[t + d];
    }
    __syncthreads();
  }
  if (!t) {
    best[blockIdx.x] = v[0];
    pos[blockIdx.x] = ix[0];
  }
}
int main(int ac, char **av) {
  int P = atoi(arg(ac, av, "--P", "3406")), n = (P + 31) / 32 * 32, L = 28,
      ns = atoi(arg(ac, av, "--samples", "128"));
  std::string qdir = arg(ac, av, "--dir", "."),
              lp = arg(ac, av, "--factor", ""),
              xp = arg(ac, av, "--inverse", "");
  H Q = readh(qdir + "/S0", P, L), HL = readh(lp, n, L);
  H HX;
  const bool inverse = !xp.empty();
  if (inverse)
    HX = readh(xp, n, L);
  if (inverse) {
    auto nq = norms(Q, P), nl = norms(HL, P), nx = norms(HX, P);
    printf(
        "NORMS log2 Q_inf %.9f L_one %.9f L_inf %.9f X_one %.9f X_inf %.9f\n",
        log2(nq.inf), log2(nl.one), log2(nl.inf), log2(nx.one), log2(nx.inf));
    printf("CONDITION_ESTIMATE log2 lower %.9f upper %.9f (stored X; not a "
           "certified bound)\n",
           log2(nq.max_diag) + log2(nx.max_col_square),
           log2(nq.inf) + log2(nx.one) + log2(nx.inf));
  }
  int perturb = atoi(arg(ac, av, "--perturb-bits", "0"));
  if (perturb) {
    int bit = 32 * L - 3 - perturb;
    if (bit < 0)
      throw std::runtime_error("perturbation too small");
    unsigned long long carry = 1ULL << (bit % 32);
    for (int t = bit / 32; t < L && carry; t++) {
      unsigned long long v = uint32_t(HL.m[size_t(t) * n * n]) + carry;
      HL.m[size_t(t) * n * n] = int32_t(v);
      carry = v >> 32;
    }
    printf("PERTURB L00 relative scale 2^-%d\n", perturb);
  }
  Gram_Options o;
  o.fused_residues = true;
  o.prime_group = 32;
  Block_Factor F(n, 1, L, 112, o);
  auto S0 = F.alloc(P, P), S = F.alloc(n, n), Lf = F.alloc(n, n);
  cudaMemcpy(S0.m, Q.m.data(), Q.m.size() * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(S0.e, Q.e.data(), Q.e.size() * 8, cudaMemcpyHostToDevice);
  mpf_zero(S, 0);
  mpf_sub_block(S, 0, 0, S0, true, 0);
  F.release(S0);
  cudaMemcpy(Lf.m, HL.m.data(), HL.m.size() * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(Lf.e, HL.e.data(), HL.e.size() * 8, cudaMemcpyHostToDevice);
  Mpf_Batch R;
  if (inverse) {
    auto X = F.alloc(n, n);
    cudaMemcpy(X.m, HX.m.data(), HX.m.size() * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(X.e, HX.e.data(), HX.e.size() * 8, cudaMemcpyHostToDevice);
    F.mm(X, false, Lf, false, R);
    F.release(X);
    // Identity is also the diagonal normalization for the inverse residual.
    std::vector<int32_t> im(size_t(L) * n * n, 0);
    std::vector<int64_t> ie(size_t(n) * n, MPF_EZ);
    for (int i = 0; i < n; i++) {
      im[size_t(L - 1) * n * n + size_t(i) * n + i] = 1 << 29;
      ie[size_t(i) * n + i] = -(32LL * L - 3);
    }
    cudaMemcpy(S.m, im.data(), im.size() * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(S.e, ie.data(), ie.size() * 8, cudaMemcpyHostToDevice);
  } else
    F.mm(Lf, false, Lf, true, R);
  mpf_batch_add(R, S, -1, R, 0);
  F.check();
  int64_t blocks = (int64_t(n) * n + 255) / 256;
  double *db;
  int64_t *di;
  cudaMalloc(&db, blocks * 8);
  cudaMalloc(&di, blocks * 8);
  reduce_res<<<blocks, 256>>>(R.m, R.e, S.m, S.e, L, n, P, db, di);
  std::vector<double> b(blocks);
  std::vector<int64_t> ix(blocks);
  cudaMemcpy(b.data(), db, blocks * 8, cudaMemcpyDeviceToHost);
  cudaMemcpy(ix.data(), di, blocks * 8, cudaMemcpyDeviceToHost);
  std::vector<int> order;
  for (int i = 0; i < blocks; i++)
    if (std::isfinite(b[i]))
      order.push_back(i);
  std::sort(order.begin(), order.end(),
            [&](int i, int j) { return b[i] > b[j]; });
  printf("%s entries %lld log2 %.9f row %lld col %lld\n",
         inverse ? "ALL_INVERSE_RESIDUAL" : "ALL_RESIDUAL",
         (long long)P * (P + 1) / 2, b[order[0]], (long long)(ix[order[0]] % n),
         (long long)(ix[order[0]] / n));
  mpf_set_default_prec(2048);
  mpf_t a, c, s, t, d;
  mpf_inits(a, c, s, t, d, nullptr);
  mpz_t z;
  mpz_init(z);
  std::mt19937 rng(325);
  double worst = -1e9, diagmin = 1e9, diagmax = -1e9;
  for (int i = 0; i < P; i++) {
    get(a, HL, i, i, z);
    diagmin = std::min(diagmin, lg(a));
    diagmax = std::max(diagmax, lg(a));
  }
  printf("DIAG_GMP min %.9f max %.9f condition_lower_bound_log2 %.9f\n",
         diagmin, diagmax, 2 * (diagmax - diagmin));
  int nw = std::min(32, int(order.size()));
  for (int r = 0; r < ns + nw; r++) {
    int i, j;
    if (r < nw) {
      i = ix[order[r]] % n;
      j = ix[order[r]] / n;
    } else {
      i = rng() % P;
      j = rng() % P;
      if (i < j)
        std::swap(i, j);
    }
    mpf_set_ui(s, 0);
    for (int k = inverse ? j : 0; k <= (inverse ? i : j); k++) {
      get(a, inverse ? HX : HL, i, k, z);
      get(c, HL, inverse ? k : j, inverse ? j : k, z);
      mpf_mul(t, a, c);
      mpf_add(s, s, t);
    }
    if (inverse)
      mpf_set_ui(a, i == j);
    else
      get(a, Q, i, j, z);
    mpf_sub(s, s, a);
    if (inverse)
      mpf_set_ui(d, 1);
    else {
      get(a, Q, i, i, z);
      get(c, Q, j, j, z);
      mpf_mul(d, a, c);
      mpf_sqrt(d, d);
    }
    double err = lg(s) - lg(d);
    worst = std::max(worst, err);
    if (r < 4)
      printf("GMP_ENTRY %d %d log2 %.9f\n", i, j, err);
  }
  printf("GMP_RESIDUAL entries %d log2 %.9f\n", ns + nw, worst);
  double limit = atof(arg(ac, av, "--max-log2", "0"));
  bool ok = b[order[0]] <= limit && worst <= limit;
  printf("CHECK %s (limit %.1f)\n", ok ? "PASSED" : "FAILED", limit);
  return ok ? 0 : 1;
}
