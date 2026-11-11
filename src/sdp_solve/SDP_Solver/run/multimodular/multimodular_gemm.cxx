#include "multimodular_blas.hxx"

#include "sdpb_util/bigint_shared_memory/fmpz/blas_residues.hxx"
#include "sdpb_util/assert.hxx"

#include <cblas.h>

#include <chrono>
#include <cstdint>
#include <limits>
#include <vector>

void Multimodular_Blas_Timings::add(const Multimodular_Blas_Timings &other)
{
  base_seconds += other.base_seconds;
  to_int_seconds += other.to_int_seconds;
  blas_seconds += other.blas_seconds;
  from_int_seconds += other.from_int_seconds;
  num_gemms += other.num_gemms;
  max_num_primes = std::max(max_num_primes, other.max_num_primes);
}

// Algorithm of multimodular_gemm_sub(C, A, op(B)):
//
// 1. Fixed point: row i of A and column j of op(B) are multiplied by
//    2^(bits - e), e = binary exponent of their largest entry, and truncated
//    to integers a_il, b_lj with |a|, |b| < 2^bits.
// 2. Residues modulo primes p_1..p_n (p < 2^21, p^2 k < 2^53, product M of
//    the primes > 8 max|c_ij|), computed by BLAS
//    (sdpb_util/bigint_shared_memory/fmpz/blas_residues.hxx). The residues
//    of A are multiplied by
//    u_p = (M/p)^-1 mod p for the CRT.
// 3. For each prime, cblas_dgemm of the residue matrices (exact in double).
// 4. CRT, also by BLAS (blas_residues.hxx). Its truncation error is
//    < 2^(bits - 89), far below the truncation error of step 1 (~k 2^bits).
// 5. C_ij -= c_ij 2^(e_i + e_j - 2 bits).

namespace
{
  using Clock = std::chrono::steady_clock;
  double seconds_since(const Clock::time_point &start)
  { return std::chrono::duration<double>(Clock::now() - start).count(); }

  using Blas_Residues::Context;
  using Blas_Residues::Crt_Workspace;
  using Blas_Residues::get_context;
  using Blas_Residues::mod_p;

  // Binary exponent e of x, such that 2^(e-1) <= |x| < 2^e.
  // Returns false for x == 0.
  bool binary_exponent(const El::BigFloat &x, long &e)
  {
    if(mpf_sgn(x.gmp_float.get_mpf_t()) == 0)
      return false;
    mpf_get_d_2exp(&e, x.gmp_float.get_mpf_t());
    return true;
  }

  // One prime's product: C = (A B mod p) in [0, p) as int32, from int32
  // residue matrices A (m x k, ld lda) and B (k x n, ld ldb), converted to
  // double for cblas_dgemm. If quotients != nullptr, adds C/p to it
  // (for the q of the BLAS CRT).
  struct Prime_Gemm_Workspace
  { std::vector<double> A, B, C; };
  void prime_gemm(const Context &context, size_t prime, int m, int n, int k,
                  const int32_t *A, size_t lda, const int32_t *B, size_t ldb,
                  int32_t *C, double *quotients, Prime_Gemm_Workspace &work)
  {
    work.A.resize(size_t(m) * k);
    work.B.resize(size_t(k) * n);
    work.C.resize(size_t(m) * n);
    for(int l = 0; l < k; ++l)
      {
        const int32_t *__restrict a = A + l * lda;
        double *__restrict ad = work.A.data() + size_t(l) * m;
        for(int i = 0; i < m; ++i)
          ad[i] = a[i];
      }
    for(int j = 0; j < n; ++j)
      {
        const int32_t *__restrict b = B + j * ldb;
        double *__restrict bd = work.B.data() + size_t(j) * k;
        for(int l = 0; l < k; ++l)
          bd[l] = b[l];
      }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0,
                work.A.data(), m, work.B.data(), k, 0.0, work.C.data(), m);
    const double p = context.primes[prime];
    const double inv_p = context.inverse_primes[prime];
    const size_t size = size_t(m) * n;
    const double *__restrict c = work.C.data();
    int32_t *__restrict out = C;
    if(quotients != nullptr)
      {
        double *__restrict q = quotients;
        for(size_t e = 0; e < size; ++e)
          {
            const double r = mod_p(c[e], p, inv_p);
            out[e] = int32_t(r);
            q[e] += r * inv_p;
          }
      }
    else
      {
        for(size_t e = 0; e < size; ++e)
          out[e] = int32_t(mod_p(c[e], p, inv_p));
      }
  }
}

namespace
{
  constexpr long no_exponent = std::numeric_limits<long>::min();

  // C_residues: num_primes blocks (stride m*width) of the m x width products
  // reduced to [0, p) (times u_p unless flint_crt); quotients[e] =
  // sum_p y_p / p (BLAS CRT only).
  // C(i, j) -= c_ij 2^(row_exp[i] + col_exp[j] - shift_base)
  // for row_start[j] <= i < m.
  void crt_subtract(const Context &context, bool flint_crt,
                    const int32_t *C_residues, const double *quotients, int m,
                    int width, El::Matrix<El::BigFloat> &C,
                    const std::vector<long> &row_exp,
                    const std::vector<long> &col_exp, long shift_base,
                    const std::vector<int> &row_start, Crt_Workspace &work,
                    const int *row_end = nullptr)
  {
    const size_t C_stride = size_t(m) * width;
    mpf_t value;
    mpf_init2(value, El::gmp::Precision());
    constexpr size_t batch = 1024;
    for(size_t e0 = 0; e0 < C_stride; e0 += batch)
      {
        const size_t nb = std::min(batch, C_stride - e0);
        if(!flint_crt)
          Blas_Residues::blas_crt_high_chunks(context, C_residues + e0,
                                              C_stride, nb, work);
        for(size_t eb = 0; eb < nb; ++eb)
          {
            const size_t e = e0 + eb;
            const int i = e % m, j = e / m;
            if(col_exp[j] == no_exponent || row_exp[i] == no_exponent
               || i < row_start[j] || (row_end != nullptr && i > row_end[j]))
              continue;
            const long shift = row_exp[i] + col_exp[j] - shift_base;
            const bool nonzero
              = flint_crt
                  ? Blas_Residues::flint_crt_to_mpf(
                      context, C_residues + e, C_stride, shift, work, value)
                  : Blas_Residues::blas_crt_to_mpf(context, work.S.data() + eb,
                                                   nb, quotients[e], shift,
                                                   work, value);
            if(!nonzero)
              continue;
            auto &c = C.Ref(i, j).gmp_float;
            mpf_sub(c.get_mpf_t(), c.get_mpf_t(), value);
          }
      }
    mpf_clear(value);
  }

  // Exponents of the largest entry of each row of A (rows = true) or of
  // each column (rows = false); no_exponent for zero rows/columns.
  template <class Get>
  std::vector<long> max_exponents(int count, int length, const Get &get)
  {
    std::vector<long> result(count, no_exponent);
    for(int a = 0; a < count; ++a)
      for(int l = 0; l < length; ++l)
        {
          long e;
          if(binary_exponent(get(a, l), e))
            result[a] = std::max(result[a], e);
        }
    return result;
  }
}

void multimodular_gemm_sub(El::Matrix<El::BigFloat> &C,
                           const El::Matrix<El::BigFloat> &A,
                           const El::Orientation orientation_B,
                           const El::Matrix<El::BigFloat> &B,
                           const bool lower_only,
                           const Multimodular_Blas_Config &config,
                           Multimodular_Blas_Timings &timings)
{
  std::vector<int> row_start(C.Width(), 0);
  if(lower_only)
    for(int j = 0; j < C.Width(); ++j)
      row_start[j] = j;
  multimodular_gemm_sub(C, A, orientation_B, B, row_start, config, timings);
}

void multimodular_gemm_sub(El::Matrix<El::BigFloat> &C,
                           const El::Matrix<El::BigFloat> &A,
                           const El::Orientation orientation_B,
                           const El::Matrix<El::BigFloat> &B,
                           const std::vector<int> &row_start,
                           const Multimodular_Blas_Config &config,
                           Multimodular_Blas_Timings &timings)
{
  multimodular_gemm_sub(C, A, orientation_B, B, row_start,
                        std::vector<int>(C.Width(), C.Height() - 1), config,
                        timings);
}

// C(i, j) -= (A op(B))(i, j) for row_start[j] <= i <= row_end[j] only.
// Each column chunk multiplies only the rows its columns need
// (the upper triangle of a distributed trailing update).
void multimodular_gemm_sub(El::Matrix<El::BigFloat> &C,
                           const El::Matrix<El::BigFloat> &A,
                           const El::Orientation orientation_B,
                           const El::Matrix<El::BigFloat> &B,
                           const std::vector<int> &row_start,
                           const std::vector<int> &row_end,
                           const Multimodular_Blas_Config &config,
                           Multimodular_Blas_Timings &timings)
{
  const bool transpose_B = orientation_B != El::NORMAL;
  const int m = A.Height();
  const int k = A.Width();
  const int n = transpose_B ? B.Height() : B.Width();
  ASSERT_EQUAL(C.Height(), m);
  ASSERT_EQUAL(C.Width(), n);
  ASSERT_EQUAL(transpose_B ? B.Width() : B.Height(), k);
  ASSERT_EQUAL(row_start.size(), size_t(n));
  ASSERT_EQUAL(row_end.size(), size_t(n));
  if(m == 0 || n == 0 || k == 0)
    return;

  auto op_B = [&](int l, int j) -> const El::BigFloat & {
    return transpose_B ? B.CRef(j, l) : B.CRef(l, j);
  };

  const long bits = El::gmp::Precision() + config.guard_bits;

  // Memory bound on the residues of A: split the rows of C and A
  if(config.gemm_A_bytes > 0 && m > 16
     && sizeof(int32_t) * size_t(m) * k * get_context(bits, bits, k).num_primes
          > config.gemm_A_bytes)
    {
      const int m1 = m / 2;
      auto C1 = C(El::IR(0, m1), El::ALL);
      auto C2 = C(El::IR(m1, m), El::ALL);
      std::vector<int> row_start2(n), row_end1(n), row_end2(n);
      for(int j = 0; j < n; ++j)
        {
          row_start2[j] = std::max(0, row_start[j] - m1);
          row_end1[j] = std::min(m1 - 1, row_end[j]);
          row_end2[j] = row_end[j] - m1;
        }
      multimodular_gemm_sub(C1, A(El::IR(0, m1), El::ALL), orientation_B, B,
                            row_start, row_end1, config, timings);
      multimodular_gemm_sub(C2, A(El::IR(m1, m), El::ALL), orientation_B, B,
                            row_start2, row_end2, config, timings);
      return;
    }

  auto start = Clock::now();

  const std::vector<long> row_exp = max_exponents(
    m, k, [&](int i, int l) -> const El::BigFloat & { return A.CRef(i, l); });
  const std::vector<long> col_exp = max_exponents(
    n, k, [&](int j, int l) -> const El::BigFloat & { return op_B(l, j); });

  const Context &context = get_context(bits, bits, k);
  const size_t num_primes = context.num_primes;
  timings.max_num_primes = std::max(timings.max_num_primes, num_primes);
  ++timings.num_gemms;
  const bool flint_crt = config.flint_crt;

  // Residues of A (times u_p for the BLAS CRT):
  // for prime p, an m x k column-major matrix at p*m*k
  const size_t A_stride = size_t(m) * k;
  std::vector<int32_t> A_residues(num_primes * A_stride);
  Blas_Residues::compute_blas_residues(
    context, A_stride,
    [&](size_t e, long &shift) -> const __mpf_struct * {
      const int i = e % m, l = e / m;
      if(row_exp[i] == no_exponent)
        return nullptr;
      shift = bits - row_exp[i];
      return A.CRef(i, l).gmp_float.get_mpf_t();
    },
    !flint_crt, A_residues.data());
  timings.to_int_seconds += seconds_since(start);

  const int chunk = std::max(1, std::min(config.column_chunk, n));
  std::vector<int32_t> B_residues(num_primes * size_t(k) * chunk);
  std::vector<int32_t> C_residues(num_primes * size_t(m) * chunk);
  std::vector<double> quotients;
  Crt_Workspace work(context);
  Prime_Gemm_Workspace gemm_work;

  for(int j0 = 0; j0 < n; j0 += chunk)
    {
      const int width = std::min(chunk, n - j0);
      const size_t B_stride = size_t(k) * width;

      start = Clock::now();
      Blas_Residues::compute_blas_residues(
        context, B_stride,
        [&](size_t e, long &shift) -> const __mpf_struct * {
          const int l = e % k, j = j0 + e / k;
          if(col_exp[j] == no_exponent)
            return nullptr;
          shift = bits - col_exp[j];
          return op_B(l, j).gmp_float.get_mpf_t();
        },
        false, B_residues.data());
      timings.to_int_seconds += seconds_since(start);

      // rows r0..r1 are the only ones this chunk's columns need
      int r0 = m, r1 = -1;
      for(int j = j0; j < j0 + width; ++j)
        if(row_start[j] <= row_end[j])
          {
            r0 = std::min(r0, std::max(0, row_start[j]));
            r1 = std::max(r1, std::min(m - 1, row_end[j]));
          }
      if(r1 < r0)
        continue;
      const int mc = r1 - r0 + 1;
      const size_t Cc_stride = size_t(mc) * width;

      start = Clock::now();
      quotients.assign(flint_crt ? 0 : Cc_stride, 0.0);
      for(size_t p = 0; p < num_primes; ++p)
        prime_gemm(context, p, mc, width, k,
                   A_residues.data() + p * A_stride + r0, m,
                   B_residues.data() + p * B_stride, k,
                   C_residues.data() + p * Cc_stride,
                   flint_crt ? nullptr : quotients.data(), gemm_work);
      timings.blas_seconds += seconds_since(start);

      start = Clock::now();
      auto C_chunk = C(El::IR(r0, r1 + 1), El::IR(j0, j0 + width));
      const std::vector<long> chunk_col_exp(col_exp.begin() + j0,
                                            col_exp.begin() + j0 + width);
      const std::vector<long> chunk_row_exp(row_exp.begin() + r0,
                                            row_exp.begin() + r1 + 1);
      std::vector<int> chunk_row_start(width), chunk_row_end(width);
      for(int j = 0; j < width; ++j)
        {
          chunk_row_start[j] = row_start[j0 + j] - r0;
          chunk_row_end[j] = row_end[j0 + j] - r0;
        }
      crt_subtract(context, flint_crt, C_residues.data(), quotients.data(),
                   mc, width, C_chunk, chunk_row_exp, chunk_col_exp,
                   2 * bits, chunk_row_start, work, chunk_row_end.data());
      timings.from_int_seconds += seconds_since(start);
    }
}

// Left-looking blocked Trsm with each conversion done once.
// For each column chunk, for each panel of h <= base_size rows:
//   B_p -= L_{p,<p} X_{<p}   (multimodular GEMM, inner dimension r0)
//   X_p = L_pp^{-1} B_p      (substitution)
//   residues of X_p are stored for the following panels.
// X is stored with one fixed-point scale per column, 2^(bits_X - e_j),
// where e_j = (largest exponent seen in column j) + growth_bits; if a
// later panel exceeds 2^e_j, e_j is raised and the stored residues of that
// column are recomputed from the BigFloat values.
void multimodular_trsm_left_looking(const El::Matrix<El::BigFloat> &L,
                                    El::Matrix<El::BigFloat> &B,
                                    const Multimodular_Blas_Config &config,
                                    Multimodular_Blas_Timings &timings)
{
  const int P = L.Height();
  const int N = B.Width();
  ASSERT_EQUAL(L.Width(), P);
  ASSERT_EQUAL(B.Height(), P);
  const int b = std::max(1, config.base_size);
  if(P == 0 || N == 0)
    return;
  if(P <= b)
    {
      auto start = Clock::now();
      substitution_trsm(L, B);
      timings.base_seconds += seconds_since(start);
      return;
    }
  const long bits_L = El::gmp::Precision() + config.guard_bits;
  const long bits_X = bits_L + config.growth_bits;
  const Context &context = get_context(bits_L, bits_X, P);
  const size_t num_primes = context.num_primes;
  timings.max_num_primes = std::max(timings.max_num_primes, num_primes);
  const bool flint_crt = config.flint_crt;

  const int chunk = std::max<size_t>(
    1, std::min<size_t>(N, config.left_looking_bytes
                             / (sizeof(double) * num_primes * P)));
  const size_t X_stride = size_t(P) * chunk; // per prime
  std::vector<int32_t> X_residues(num_primes * X_stride);
  std::vector<int32_t> A_residues, C_residues, panel_residues;
  std::vector<double> quotients;
  Crt_Workspace work(context);
  Prime_Gemm_Workspace gemm_work;

  for(int j0 = 0; j0 < N; j0 += chunk)
    {
      const int width = std::min(chunk, N - j0);
      std::vector<long> col_e(width, no_exponent);
      for(int r0 = 0; r0 < P; r0 += b)
        {
          const int h = std::min(b, P - r0);
          auto R = B(El::IR(r0, r0 + h), El::IR(j0, j0 + width));
          if(r0 > 0)
            {
              auto start = Clock::now();
              const auto A = L(El::IR(r0, r0 + h), El::IR(0, r0));
              const std::vector<long> row_exp = max_exponents(
                h, r0, [&](int i, int l) -> const El::BigFloat & {
                  return A.CRef(i, l);
                });
              const size_t A_stride = size_t(h) * r0;
              A_residues.resize(num_primes * A_stride);
              Blas_Residues::compute_blas_residues(
                context, A_stride,
                [&](size_t e, long &shift) -> const __mpf_struct * {
                  const int i = e % h, l = e / h;
                  if(row_exp[i] == no_exponent)
                    return nullptr;
                  shift = bits_L - row_exp[i];
                  return A.CRef(i, l).gmp_float.get_mpf_t();
                },
                !flint_crt, A_residues.data());
              timings.to_int_seconds += seconds_since(start);

              start = Clock::now();
              const size_t C_stride = size_t(h) * width;
              C_residues.resize(num_primes * C_stride);
              quotients.assign(flint_crt ? 0 : C_stride, 0.0);
              for(size_t p = 0; p < num_primes; ++p)
                prime_gemm(context, p, h, width, r0,
                           A_residues.data() + p * A_stride, h,
                           X_residues.data() + p * X_stride, P,
                           C_residues.data() + p * C_stride,
                           flint_crt ? nullptr : quotients.data(), gemm_work);
              timings.blas_seconds += seconds_since(start);
              ++timings.num_gemms;

              start = Clock::now();
              crt_subtract(context, flint_crt, C_residues.data(),
                           quotients.data(), h, width, R, row_exp, col_e,
                           bits_L + bits_X, std::vector<int>(width, 0), work);
              timings.from_int_seconds += seconds_since(start);
            }

          auto start = Clock::now();
          substitution_trsm(L(El::IR(r0, r0 + h), El::IR(r0, r0 + h)), R);
          timings.base_seconds += seconds_since(start);

          start = Clock::now();
          // Column scales: first nonzero panel sets them; growth beyond the
          // headroom raises them and recomputes the stored residues.
          const std::vector<long> panel_exp = max_exponents(
            width, h, [&](int j, int i) -> const El::BigFloat & {
              return R.CRef(i, j);
            });
          for(int j = 0; j < width; ++j)
            {
              if(panel_exp[j] == no_exponent || panel_exp[j] <= col_e[j])
                continue;
              const bool had_nonzero = col_e[j] != no_exponent;
              col_e[j] = panel_exp[j] + config.growth_bits;
              if(had_nonzero && r0 > 0)
                {
                  ++timings.num_rescales;
                  Blas_Residues::compute_blas_residues(
                    context, r0,
                    [&](size_t e, long &shift) -> const __mpf_struct * {
                      shift = bits_X - col_e[j];
                      return B.CRef(e, j0 + j).gmp_float.get_mpf_t();
                    },
                    false, X_residues.data() + size_t(P) * j, X_stride);
                }
            }
          // Residues of the new panel rows, then into the X store
          const size_t panel_stride = size_t(h) * width;
          panel_residues.resize(num_primes * panel_stride);
          Blas_Residues::compute_blas_residues(
            context, panel_stride,
            [&](size_t e, long &shift) -> const __mpf_struct * {
              const int i = e % h, j = e / h;
              if(col_e[j] == no_exponent)
                return nullptr;
              shift = bits_X - col_e[j];
              return R.CRef(i, j).gmp_float.get_mpf_t();
            },
            false, panel_residues.data());
          for(size_t p = 0; p < num_primes; ++p)
            for(int j = 0; j < width; ++j)
              std::copy_n(
                panel_residues.data() + p * panel_stride + size_t(j) * h, h,
                X_residues.data() + p * X_stride + size_t(j) * P + r0);
          timings.to_int_seconds += seconds_since(start);
        }
    }
}

void check_multimodular_precision(const Multimodular_Blas_Config &config,
                                  const long precision)
{
  // As in multimodular_gemm() and the left-looking Trsm below
  const long bits = precision + config.guard_bits;
  const long bits_X = bits + config.growth_bits;
  for(slong k = 1; k <= (slong(1) << 24); k *= 2)
    if(!Blas_Residues::is_exact(bits, bits, k)
       || !Blas_Residues::is_exact(bits, bits_X, k))
      RUNTIME_ERROR("--multimodularTrsm does not support --precision=",
                    precision,
                    ": the residues and CRT in double would not be exact "
                    "(the limit is 3200 bits).");
}
