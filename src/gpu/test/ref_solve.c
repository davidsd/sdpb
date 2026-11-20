// ref_solve: an independent high-precision
// reference for P = L^-1 B on the real Schur blocks, L = chol(S).
//
// MPFR (not GMP mpf, which the backend's own checks use), written from SDPviBe's
// documented MPF layout (sdpvibe/mpf.py: m int32 (L, *shape) limb planes, little
// endian, two's complement; e int64 (*shape); value = m * 2^e; zero: e = EZ).
// Left-looking Cholesky on the packed lower triangle and forward substitution,
// every operation rounded to `prec` bits.  For each selected column the errors
// of the backend's P and of SDPviBe's P are max_i |P_i - y_i| / max_i |y_i|.
//
//   ref_solve <S dir> <chunk> <cols.bin> <prec> [--blocks b0,b1,..] [--tri lower|upper]
//   ref_solve --selftest n prec
// cols.bin (written by q_check): "QCOL", int32 nb, P, N, L, K; int32 col[nb*K],
// kind[nb*K]; then per (b, k): B, Pb, Psv columns, P entries each, every entry
// L int32 limbs then an int64 exponent.
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>

#define EZ (-(1LL << 40))

// one block wanted: parallel inside its Cholesky instead of across blocks
static int g_inner = 0;
static int g_failed = 0;

static void mpf_entry(mpfr_t x, const uint32_t *limbs, int L, int64_t e, mpz_t z)
{
  if(e == EZ)
    {
      mpfr_set_zero(x, 1);
      return;
    }
  mpz_import(z, L, -1, 4, 0, 0, limbs);
  if(limbs[L - 1] >> 31) // two's complement: subtract 2^(32 L)
    {
      mpz_t p;
      mpz_init(p);
      mpz_setbit(p, 32 * (mp_bitcnt_t)L);
      mpz_sub(z, z, p);
      mpz_clear(p);
    }
  int r = mpfr_set_z_2exp(x, z, e, MPFR_RNDN);
  if(r != 0)
    {
      fprintf(stderr, "inexact conversion (precision too low)\n");
      exit(3);
    }
}

static double lg2(const mpfr_t x)
{
  if(mpfr_zero_p(x))
    return -1e9;
  long ex;
  double d = mpfr_get_d_2exp(&ex, x, MPFR_RNDN);
  return log2(fabs(d)) + (double)ex;
}

typedef struct
{
  int P;
  mpfr_t *Lr; // packed rows: row i at Lr[i(i+1)/2 .. +i]
} chol_t;

#define LIDX(i, k) ((int64_t)(i) * ((i) + 1) / 2 + (k))

// S given as a callback filling the packed lower triangle; returns the
// smallest pivot's log2 (and -1e9 if a pivot is not positive)
static double chol(chol_t *C, mpfr_prec_t prec)
{
  const int P = C->P;
  mpfr_t s, t;
  mpfr_inits2(prec, s, t, (mpfr_ptr)0);
  double pmin = 1e9;
  for(int j = 0; j < P; ++j)
    {
      mpfr_t *Lj = C->Lr + LIDX(j, 0);
      // pivot
      mpfr_set(s, Lj[j], MPFR_RNDN);
      for(int k = 0; k < j; ++k)
        {
          mpfr_sqr(t, Lj[k], MPFR_RNDN);
          mpfr_sub(s, s, t, MPFR_RNDN);
        }
      if(mpfr_sgn(s) <= 0)
        {
          fprintf(stderr, "pivot %d not positive\n", j);
          pmin = -1e9;
          mpfr_abs(s, s, MPFR_RNDN);
        }
      else
        pmin = fmin(pmin, lg2(s));
      mpfr_sqrt(Lj[j], s, MPFR_RNDN);
      // column j below the pivot: L_ij = (S_ij - sum_k L_ik L_jk) / L_jj
#pragma omp parallel for schedule(dynamic, 8) if(g_inner)
      for(int i = j + 1; i < P; ++i)
        {
          mpfr_t *Li = C->Lr + LIDX(i, 0);
          mpfr_t a, b;
          mpfr_inits2(prec, a, b, (mpfr_ptr)0);
          mpfr_set(a, Li[j], MPFR_RNDN);
          for(int k = 0; k < j; ++k)
            {
              mpfr_mul(b, Li[k], Lj[k], MPFR_RNDN);
              mpfr_sub(a, a, b, MPFR_RNDN);
            }
          mpfr_div(Li[j], a, Lj[j], MPFR_RNDN);
          mpfr_clears(a, b, (mpfr_ptr)0);
        }
    }
  mpfr_clears(s, t, (mpfr_ptr)0);
  return pmin;
}

// y = L^-1 b (b overwritten by y)
static void forward_sub(const chol_t *C, mpfr_t *y, mpfr_prec_t prec)
{
  mpfr_t t;
  mpfr_init2(t, prec);
  for(int i = 0; i < C->P; ++i)
    {
      const mpfr_t *Li = (const mpfr_t *)(C->Lr + LIDX(i, 0));
      for(int k = 0; k < i; ++k)
        {
          mpfr_mul(t, Li[k], y[k], MPFR_RNDN);
          mpfr_sub(y[i], y[i], t, MPFR_RNDN);
        }
      mpfr_div(y[i], y[i], Li[i], MPFR_RNDN);
    }
  mpfr_clear(t);
}

// relative error of computed column c against y: max|c - y| / max|y| (log2)
static double colerr(mpfr_t *c, mpfr_t *y, int P, mpfr_prec_t prec)
{
  mpfr_t d, my, md;
  mpfr_inits2(prec, d, my, md, (mpfr_ptr)0);
  mpfr_set_zero(my, 1);
  mpfr_set_zero(md, 1);
  for(int i = 0; i < P; ++i)
    {
      mpfr_sub(d, c[i], y[i], MPFR_RNDN);
      mpfr_abs(d, d, MPFR_RNDN);
      if(mpfr_cmp(d, md) > 0)
        mpfr_set(md, d, MPFR_RNDN);
      mpfr_abs(d, y[i], MPFR_RNDN);
      if(mpfr_cmp(d, my) > 0)
        mpfr_set(my, d, MPFR_RNDN);
    }
  double r = lg2(md) - lg2(my);
  if(mpfr_zero_p(md))
    r = -1e9;
  mpfr_clears(d, my, md, (mpfr_ptr)0);
  return r;
}

static double colmax(mpfr_t *c, int P)
{
  double m = -1e9;
  for(int i = 0; i < P; ++i)
    m = fmax(m, lg2(c[i]));
  return m;
}

static void *map(const char *path, size_t *len)
{
  int fd = open(path, O_RDONLY);
  if(fd < 0)
    {
      perror(path);
      exit(2);
    }
  struct stat st;
  fstat(fd, &st);
  *len = st.st_size;
  void *p = mmap(0, *len, PROT_READ, MAP_SHARED, fd, 0);
  if(p == MAP_FAILED)
    {
      perror("mmap");
      exit(2);
    }
  close(fd);
  return p;
}

// ---- self test: S = L0 L0^T and B = L0 Y with small integers, so P = Y exactly
static int selftest(int n, mpfr_prec_t prec)
{
  chol_t C = {n, malloc(sizeof(mpfr_t) * LIDX(n, 0))};
  srand(7);
  long *L0 = calloc((size_t)n * n, sizeof(long)), *Y = malloc(sizeof(long) * n);
  for(int i = 0; i < n; ++i)
    for(int k = 0; k <= i; ++k)
      L0[(size_t)i * n + k] = (k == i) ? 1 + rand() % 50 : (rand() % 2001) - 1000;
  for(int i = 0; i < n; ++i)
    Y[i] = (rand() % 2001) - 1000;
  mpz_t acc, t;
  mpz_inits(acc, t, NULL);
  for(int i = 0; i < n; ++i)
    for(int j = 0; j <= i; ++j)
      {
        mpz_set_ui(acc, 0);
        for(int k = 0; k <= j; ++k)
          {
            mpz_set_si(t, L0[(size_t)i * n + k]);
            mpz_mul_si(t, t, L0[(size_t)j * n + k]);
            mpz_add(acc, acc, t);
          }
        mpfr_init2(C.Lr[LIDX(i, j)], prec);
        if(mpfr_set_z(C.Lr[LIDX(i, j)], acc, MPFR_RNDN))
          return printf("selftest: S not exact\n"), 1;
      }
  mpfr_t *y = malloc(sizeof(mpfr_t) * n);
  for(int i = 0; i < n; ++i)
    {
      mpz_set_ui(acc, 0);
      for(int k = 0; k <= i; ++k)
        {
          mpz_set_si(t, L0[(size_t)i * n + k]);
          mpz_mul_si(t, t, Y[k]);
          mpz_add(acc, acc, t);
        }
      mpfr_init2(y[i], prec);
      mpfr_set_z(y[i], acc, MPFR_RNDN);
    }
  chol(&C, prec);
  // L = L0 D^-1/2 with... no: chol(L0 L0^T) = L0 exactly (positive diagonal), so y = Y
  forward_sub(&C, y, prec);
  mpfr_t *yy = malloc(sizeof(mpfr_t) * n);
  for(int i = 0; i < n; ++i)
    {
      mpfr_init2(yy[i], prec);
      mpfr_set_si(yy[i], Y[i], MPFR_RNDN);
    }
  double e = colerr(y, yy, n, prec), eL = -1e9;
  mpfr_t d;
  mpfr_init2(d, prec);
  for(int i = 0; i < n; ++i)
    for(int k = 0; k <= i; ++k)
      {
        mpfr_sub_si(d, C.Lr[LIDX(i, k)], L0[(size_t)i * n + k], MPFR_RNDN);
        eL = fmax(eL, lg2(d));
      }
  // reader: -3 * 2^-5 and (2^893 + 1) * 2^-900 in 28 two's complement limbs
  {
    uint32_t w[28];
    mpz_t zz;
    mpz_init(zz);
    mpfr_t v, r;
    mpfr_inits2(prec, v, r, (mpfr_ptr)0);
    for(int q = 0; q < 28; ++q)
      w[q] = 0xffffffffu;
    w[0] = 0xfffffffdu;
    mpf_entry(v, w, 28, -5, zz);
    mpfr_set_si_2exp(r, -3, -5, MPFR_RNDN);
    int ok1 = mpfr_equal_p(v, r);
    memset(w, 0, sizeof w);
    w[0] = 1;
    w[27] = 1u << 29;
    mpf_entry(v, w, 28, -900, zz);
    mpfr_set_ui_2exp(r, 1, 893, MPFR_RNDN);
    mpfr_add_ui(r, r, 1, MPFR_RNDN);
    mpfr_mul_2si(r, r, -900, MPFR_RNDN);
    int ok2 = mpfr_equal_p(v, r);
    mpf_entry(v, w, 28, EZ, zz);
    printf("selftest reader: %s %s %s\n", ok1 ? "ok" : "FAIL", ok2 ? "ok" : "FAIL", mpfr_zero_p(v) ? "ok" : "FAIL");
  }
  printf("selftest n=%d prec=%ld: max|L - L0| = 2^%.1f, max|y - Y|/max|Y| = 2^%.1f\n", n, (long)prec, eL, e);
  return 0;
}

int main(int argc, char **argv)
{
  if(argc >= 4 && !strcmp(argv[1], "--selftest"))
    return selftest(atoi(argv[2]), atol(argv[3]));
  if(argc < 5)
    {
      fprintf(stderr, "usage: ref_solve <S dir> <chunk> <cols.bin> <prec> [--blocks b,..] [--tri lower|upper]\n");
      return 1;
    }
  const char *dir = argv[1], *chunk = argv[2], *colsf = argv[3];
  const mpfr_prec_t prec = atol(argv[4]);
  const char *blist = NULL, *save_ref = NULL, *load_ref = NULL, *extra_ref = NULL;
  int upper = 0, require_no_worse = 0;
  double max_error = INFINITY;
  // negative control: S[i][j] (and S[j][i]) of block pb multiplied by 1 + 2^-pbits
  int pb = -1, pi = 0, pj = 0, pbits = 0;
  for(int a = 5; a + 1 < argc; ++a)
    {
      if(!strcmp(argv[a], "--require-no-worse-sv")) require_no_worse=atoi(argv[a+1]);
      if(!strcmp(argv[a], "--max-log2")) max_error=atof(argv[a+1]);
      if(!strcmp(argv[a], "--extra-ref")) extra_ref = argv[a + 1];
      if(!strcmp(argv[a], "--save-ref")) save_ref = argv[a + 1];
      if(!strcmp(argv[a], "--load-ref")) load_ref = argv[a + 1];
      if(!strcmp(argv[a], "--perturb"))
        sscanf(argv[a + 1], "%d,%d,%d,%d", &pb, &pi, &pj, &pbits);
      if(!strcmp(argv[a], "--blocks"))
        blist = argv[a + 1];
      if(!strcmp(argv[a], "--tri"))
        upper = !strcmp(argv[a + 1], "upper");
    }
  if(load_ref && (save_ref || pb>=0)) {
    fprintf(stderr,"load-ref cannot be combined with save-ref or perturb\n");return 2;
  }
  size_t lc, lm, le;
  const char *cb = map(colsf, &lc);
  if(memcmp(cb, "QCOL", 4))
    return fprintf(stderr, "bad cols file\n"), 2;
  const int32_t *h = (const int32_t *)(cb + 4);
  const int nb = h[0], P = h[1], N = h[2], L = h[3], K = h[4];
  const int32_t *col = h + 5, *kind = col + (size_t)nb * K;
  const char *ent = (const char *)(kind + (size_t)nb * K);
  const size_t esz = 4 * (size_t)L + 8;
  char path[4096];
  snprintf(path, sizeof path, "%s/S%s.m.bin", dir, chunk);
  const uint32_t *Sm = map(path, &lm);
  snprintf(path, sizeof path, "%s/S%s.e.bin", dir, chunk);
  const int64_t *Se = map(path, &le);
  const int64_t plane = (int64_t)nb * P * P;
  if(lm != (size_t)L * plane * 4 || le != (size_t)plane * 8)
    return fprintf(stderr, "S size mismatch: %zu %zu vs nb=%d P=%d L=%d\n", lm, le, nb, P, L), 2;
  printf("ref_solve: %s chunk %s, nb=%d P=%d N=%d L=%d, %d columns per block, prec %ld, %s triangle\n",
         dir, chunk, nb, P, N, L, K, (long)prec, upper ? "upper" : "lower");
  int want[4096] = {0};
  if(blist)
    {
      for(const char *q = blist; *q;)
        {
          want[atoi(q)] = 1;
          while(*q && *q != ',')
            ++q;
          if(*q == ',')
            ++q;
        }
    }
  else
    for(int b = 0; b < nb; ++b)
      want[b] = 1;
  {
    int nw = 0;
    for(int b = 0; b < nb; ++b)
      nw += want[b];
    g_inner = nw == 1;
  }
  // exact symmetry of the stored blocks (limbs and exponents), all blocks
  long asym = 0;
  for(int b = 0; b < nb; ++b)
    for(int i = 0; i < P; ++i)
      for(int j = 0; j < i; ++j)
        {
          int64_t p = (int64_t)b * P * P + (int64_t)i * P + j, q = (int64_t)b * P * P + (int64_t)j * P + i;
          int same = Se[p] == Se[q];
          for(int t = 0; same && t < L; ++t)
            same = Sm[t * plane + p] == Sm[t * plane + q];
          asym += !same;
        }
  printf("S stored exactly symmetric: %s (%ld asymmetric pairs)\n", asym ? "NO" : "yes", asym);
  fflush(stdout);
#pragma omp parallel for schedule(dynamic, 1) if(!g_inner)
  for(int b = 0; b < nb; ++b)
    {
      if(!want[b])
        continue;
      double t0 = omp_get_wtime();
      chol_t C = {P, load_ref ? NULL : malloc(sizeof(mpfr_t) * LIDX(P, 0))};
      mpz_t z;
      mpz_init(z);
      uint32_t limbs[64];
      for(int i = 0; !load_ref && i < P; ++i)
        for(int j = 0; j <= i; ++j)
          {
            // lower: S[i][j] (row i, col j of SDPviBe's row-major block); upper: S[j][i]
            int64_t p = (int64_t)b * P * P + (upper ? (int64_t)j * P + i : (int64_t)i * P + j);
            for(int t = 0; t < L; ++t)
              limbs[t] = Sm[t * plane + p];
            mpfr_init2(C.Lr[LIDX(i, j)], prec);
            mpf_entry(C.Lr[LIDX(i, j)], limbs, L, Se[p], z);
          }
      if(b == pb && !load_ref)
        {
          const int ii = pi > pj ? pi : pj, jj = pi > pj ? pj : pi;
          mpfr_t t;
          mpfr_init2(t, prec);
          mpfr_mul_2si(t, C.Lr[LIDX(ii, jj)], -pbits, MPFR_RNDN);
          mpfr_add(C.Lr[LIDX(ii, jj)], C.Lr[LIDX(ii, jj)], t, MPFR_RNDN);
          mpfr_clear(t);
#pragma omp critical
          printf("NEGATIVE CONTROL: block %d, S[%d][%d] *= 1 + 2^-%d\n", b, ii, jj, pbits);
        }
      FILE *rf = NULL; mpfr_t *refy = NULL; unsigned char *seen=NULL;
      if(save_ref || load_ref) {
        if(nb != 1) { fprintf(stderr, "cached reference requires one block\n"); exit(2); }
        rf = fopen(save_ref ? save_ref : load_ref, save_ref ? "w" : "r");
        if(!rf) { perror("reference file"); exit(2); }
      }
      if(load_ref) {
        seen=calloc(N,1);
        refy=malloc(sizeof(mpfr_t)*(size_t)N*P);
        for(size_t j=0;j<(size_t)N*P;j++) mpfr_init2(refy[j],prec);
        int c,rc;
        while((rc=fscanf(rf," %d",&c))==1) {
          if(c<0 || c>=N || seen[c]) { fprintf(stderr,"bad/duplicate reference column\n"); exit(2); }
          seen[c]=1;
          for(int i=0;i<P;i++) if(!mpfr_inp_str(refy[(size_t)c*P+i],rf,16,MPFR_RNDN) || !mpfr_number_p(refy[(size_t)c*P+i])) {fprintf(stderr,"bad reference value\n"); exit(2);}
        }
        if(rc!=EOF) { fprintf(stderr,"malformed reference cache\n");exit(2); }
      }
      double pmin = load_ref ? 0 : chol(&C, prec);
      if(!load_ref && pmin <= -1e9) {
#pragma omp atomic write
        g_failed=1;
      }
      if(!load_ref) {
        double lmin = 1e9, lmax = -1e9;
        for(int i=0;i<P;i++) { double v=lg2(C.Lr[LIDX(i,i)]); lmin=fmin(lmin,v); lmax=fmax(lmax,v); }
        printf("DIAG log2 L min %.9f max %.9f condition_lower_bound_log2 %.9f\n", lmin,lmax,2*(lmax-lmin));
        fflush(stdout);
      }
      double t1 = omp_get_wtime();
      mpfr_t *y = malloc(sizeof(mpfr_t) * P), *pb = malloc(sizeof(mpfr_t) * P), *ps = malloc(sizeof(mpfr_t) * P);
      for(int i = 0; i < P; ++i)
        {
          mpfr_init2(y[i], prec);
          mpfr_init2(pb[i], prec);
          mpfr_init2(ps[i], prec);
        }
      char line[256 * 64];
      int off = 0;
      double wb = -1e9, ws = -1e9;
      for(int k = 0; k < K; ++k)
        {
          const int c = col[(size_t)b * K + k];
          if(c < 0)
            continue;
          const char *base = ent + ((size_t)b * K + k) * 3 * (size_t)P * esz;
          for(int a = 0; a < 3; ++a)
            for(int i = 0; i < P; ++i)
              {
                const char *q = base + ((size_t)a * P + i) * esz;
                memcpy(limbs, q, 4 * L);
                int64_t e;
                memcpy(&e, q + 4 * L, 8);
                mpf_entry(a == 0 ? y[i] : a == 1 ? pb[i] : ps[i], limbs, L, e, z);
              }
          double bmax = colmax(y, P);
          if(load_ref) {
            if(c>=N || !seen[c]) {fprintf(stderr,"missing reference column %d\n",c);exit(2);}
            for(int i=0;i<P;i++) mpfr_set(y[i],refy[(size_t)c*P+i],MPFR_RNDN);
          } else forward_sub(&C, y, prec);
          if(save_ref) { fprintf(rf,"%d\n",c); for(int i=0;i<P;i++) { mpfr_out_str(rf,16,0,y[i],MPFR_RNDN); fputc('\n',rf); } }
          double eb = colerr(pb, y, P, prec), es = colerr(ps, y, P, prec), ymax = colmax(y, P);
          if(!isfinite(eb) || !isfinite(es) || eb>max_error || (require_no_worse && eb>es)) {
#pragma omp atomic write
            g_failed=1;
          }
          wb = fmax(wb, eb);
          ws = fmax(ws, es);
          if(off < (int)sizeof line - 256)
            off += snprintf(line + off, sizeof line - off, "COL %d %d %d %.2f %.2f %.2f %.2f\n", b, c,
                            kind[(size_t)b * K + k], eb, es, ymax, bmax);
        }
#pragma omp critical
      {
        printf("%s", line);
        if(load_ref)
          printf("BLOCK %d: worst backend 2^%.2f, worst SDPviBe 2^%.2f; cached reference, load %.0f s, compare %.0f s\n",
                 b, wb, ws, t1 - t0, omp_get_wtime() - t1);
        else
          printf("BLOCK %d: worst backend 2^%.2f, worst SDPviBe 2^%.2f; min pivot 2^%.1f; chol %.0f s, solves %.0f s\n",
                 b, wb, ws, pmin, t1 - t0, omp_get_wtime() - t1);
        fflush(stdout);
      }
      if(extra_ref && !load_ref) {
        FILE *ex=fopen(extra_ref,"w"); if(!ex) {perror("extra reference");exit(2);}
        for(int c=0;c<128;c++) {
          int j=c<32 ? c : 32+(int64_t)(P-33)*(c-32)/95;
          for(int i=0;i<P;i++) mpfr_set_ui(y[i],i==j,MPFR_RNDN);
          forward_sub(&C,y,prec);
          fprintf(ex,"%d\n",c);
          for(int i=0;i<P;i++) {mpfr_out_str(ex,16,0,y[i],MPFR_RNDN);fputc('\n',ex);}
        }
        fclose(ex); printf("EXTRA_REFERENCE 128 deterministic inverse columns saved, total %.1f s\n",omp_get_wtime()-t0);fflush(stdout);
      }
      for(int i = 0; i < P; ++i)
        {
          mpfr_clear(y[i]);
          mpfr_clear(pb[i]);
          mpfr_clear(ps[i]);
        }
      free(y);
      free(pb);
      free(ps);
      for(int64_t q = 0; !load_ref && q < LIDX(P, 0); ++q)
        mpfr_clear(C.Lr[q]);
      free(C.Lr);
      if(rf) fclose(rf);
      if(refy) { for(size_t j=0;j<(size_t)N*P;j++) mpfr_clear(refy[j]); free(refy); }
      free(seen);
      mpz_clear(z);
    }
  if(require_no_worse || isfinite(max_error)) printf("MPFR_GATE %s\n",g_failed ? "FAILED" : "PASSED");
  return g_failed ? 1 : 0;
}
