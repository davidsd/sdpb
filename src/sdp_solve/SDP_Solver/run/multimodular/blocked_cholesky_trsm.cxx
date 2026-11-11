#include "multimodular_blas.hxx"

#include "sdpb_util/assert.hxx"

#include <cstdlib>

#include <chrono>
#include <vector>

namespace
{
  using Clock = std::chrono::steady_clock;
  double seconds_since(const Clock::time_point &start)
  { return std::chrono::duration<double>(Clock::now() - start).count(); }
}

// B := L^{-1} B by forward substitution, column by column, with
// preallocated mpf temporaries (El::Trsm on BigFloat allocates a temporary
// per multiply-add).
void substitution_trsm(const El::Matrix<El::BigFloat> &L,
                       El::Matrix<El::BigFloat> &B)
{
  const int n = L.Height();
  const mp_bitcnt_t precision = El::gmp::Precision();
  mpf_t t;
  mpf_init2(t, precision);
  // 1/L_ii: one division per row instead of one per entry
  std::vector<El::BigFloat> inverse_diagonal(n);
  for(int i = 0; i < n; ++i)
    mpf_ui_div(inverse_diagonal[i].gmp_float.get_mpf_t(), 1,
               L.CRef(i, i).gmp_float.get_mpf_t());
  for(int j = 0; j < B.Width(); ++j)
    for(int i = 0; i < n; ++i)
      {
        mpf_ptr x_i = B.Ref(i, j).gmp_float.get_mpf_t();
        for(int k = 0; k < i; ++k)
          {
            mpf_srcptr l_ik = L.CRef(i, k).gmp_float.get_mpf_t();
            mpf_srcptr x_k = B.CRef(k, j).gmp_float.get_mpf_t();
            if(mpf_sgn(l_ik) == 0 || mpf_sgn(x_k) == 0)
              continue;
            mpf_mul(t, l_ik, x_k);
            mpf_sub(x_i, x_i, t);
          }
        mpf_mul(x_i, x_i, inverse_diagonal[i].gmp_float.get_mpf_t());
      }
  mpf_clear(t);
}

// Recursive: L = [L11 0; L21 L22], B = [B1; B2]
//   B1 := L11^{-1} B1
//   B2 -= L21 B1        (multimodular GEMM)
//   B2 := L22^{-1} B2
void multimodular_trsm_left_lower(const El::Matrix<El::BigFloat> &L,
                                  El::Matrix<El::BigFloat> &B,
                                  const Multimodular_Blas_Config &config,
                                  Multimodular_Blas_Timings &timings)
{
  if(config.left_looking
     && B.Width() <= config.left_looking_max_ratio * L.Height())
    {
      multimodular_trsm_left_looking(L, B, config, timings);
      return;
    }
  const int n = L.Height();
  ASSERT_EQUAL(L.Width(), n);
  ASSERT_EQUAL(B.Height(), n);
  if(n <= config.base_size)
    {
      auto start = Clock::now();
      substitution_trsm(L, B);
      timings.base_seconds += seconds_since(start);
      return;
    }
  const int n1 = n / 2;
  const auto top = El::IR(0, n1), bottom = El::IR(n1, n), all = El::ALL;
  auto B1 = B(top, all);
  auto B2 = B(bottom, all);
  multimodular_trsm_left_lower(L(top, top), B1, config, timings);
  multimodular_gemm_sub(B2, L(bottom, top), El::NORMAL, B1, false, config,
                        timings);
  multimodular_trsm_left_lower(L(bottom, bottom), B2, config, timings);
}

// Cholesky of a small matrix in BigFloat, reading and writing only the
// lower triangle (the strictly upper part may hold anything).
void base_cholesky_lower(El::Matrix<El::BigFloat> &A)
{
  const int n = A.Height();
  const mp_bitcnt_t precision = El::gmp::Precision();
  mpf_t t, inverse;
  mpf_init2(t, precision);
  mpf_init2(inverse, precision);
  for(int j = 0; j < n; ++j)
    {
      mpf_ptr a_jj = A.Ref(j, j).gmp_float.get_mpf_t();
      for(int k = 0; k < j; ++k)
        {
          mpf_srcptr a_jk = A.CRef(j, k).gmp_float.get_mpf_t();
          mpf_mul(t, a_jk, a_jk);
          mpf_sub(a_jj, a_jj, t);
        }
      if(mpf_sgn(a_jj) <= 0)
        {
          mpf_clear(t);
          mpf_clear(inverse);
          RUNTIME_ERROR("Cholesky: matrix is not positive definite (pivot ", j,
                        " of ", n, ")");
        }
      mpf_sqrt(a_jj, a_jj);
      mpf_ui_div(inverse, 1, a_jj);
      for(int i = j + 1; i < n; ++i)
        {
          mpf_ptr a_ij = A.Ref(i, j).gmp_float.get_mpf_t();
          for(int k = 0; k < j; ++k)
            {
              mpf_mul(t, A.CRef(i, k).gmp_float.get_mpf_t(),
                      A.CRef(j, k).gmp_float.get_mpf_t());
              mpf_sub(a_ij, a_ij, t);
            }
          mpf_mul(a_ij, a_ij, inverse);
        }
    }
  mpf_clear(t);
  mpf_clear(inverse);
}

// Recursive: A = [A11 *; A21 A22]
//   A11 := L11 = chol(A11)
//   A21 := L21 = A21 L11^{-T}, i.e. L21^T = L11^{-1} A21^T
//   A22 -= L21 L21^T    (multimodular GEMM, lower triangle only)
//   A22 := chol(A22)
void multimodular_cholesky_lower(El::Matrix<El::BigFloat> &A,
                                 const Multimodular_Blas_Config &config,
                                 Multimodular_Blas_Timings &timings)
{
  const int n = A.Height();
  ASSERT_EQUAL(A.Width(), n);
  if(n <= config.base_size)
    {
      auto start = Clock::now();
      base_cholesky_lower(A);
      timings.base_seconds += seconds_since(start);
      return;
    }
  const int n1 = n / 2;
  const auto top = El::IR(0, n1), bottom = El::IR(n1, n);
  auto A11 = A(top, top);
  auto A21 = A(bottom, top);
  auto A22 = A(bottom, bottom);

  multimodular_cholesky_lower(A11, config, timings);

  El::Matrix<El::BigFloat> A21_T;
  El::Transpose(A21, A21_T);
  multimodular_trsm_left_lower(A11, A21_T, config, timings);
  for(int j = 0; j < A21.Width(); ++j)
    for(int i = 0; i < A21.Height(); ++i)
      A21.Set(i, j, A21_T.CRef(j, i));

  multimodular_gemm_sub(A22, A21, El::TRANSPOSE, A21, true, config, timings);
  multimodular_cholesky_lower(A22, config, timings);
}

void apply_multimodular_test_overrides(Multimodular_Blas_Config &config)
{
  const auto read = [](const char *name, int &value) {
    const char *text = std::getenv(name);
    if(text == nullptr)
      return;
    const int v = std::atoi(text);
    ASSERT(v > 0, name, "=", text, " must be a positive integer");
    value = v;
  };
  read("SDPB_MULTIMODULAR_BASE_SIZE", config.base_size);
  read("SDPB_MULTIMODULAR_DISTRIBUTED_MIN", config.distributed_min_size);
}
