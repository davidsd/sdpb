#include "pmp/Polynomial.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <El.hpp>
#include <vector>
#include <boost/math/tools/polynomial.hpp>
#include <mps/mps.h>

struct MPSolve_Root
{
  [[nodiscard]]
  MPSolve_Root(const El::Complex<El::BigFloat> &value,
               const El::BigFloat &radius, const mps_root_status status)
      : value(value), radius(radius), status(status)
  {}
  El::Complex<El::BigFloat> value;
  El::BigFloat radius;
  mps_root_status status;
};

std::vector<MPSolve_Root>
find_polynomial_roots(const std::vector<El::BigFloat> &polynomial_coeffs,
                      Timers &timers)
{
  if(polynomial_coeffs.size() <= 1)
    return {};

  std::vector<MPSolve_Root> result;

  const auto degree = polynomial_coeffs.size() - 1;
  const auto prec = El::gmp::Precision();

  mps_context *ctx = mps_context_new();
  mps_context_set_input_prec(ctx, prec);
  mps_context_set_output_prec(ctx, prec);

  mps_monomial_poly *mps_poly = mps_monomial_poly_new(ctx, degree);
  {
    El::BigFloat zero(0);
    mpc_t coeff;
    mpc_init2(coeff, prec);
    for(int i = 0; i < polynomial_coeffs.size(); i++)
      {
        mpf_set(mpc_Re(coeff), polynomial_coeffs.at(i).gmp_float.get_mpf_t());
        mpf_set(mpc_Im(coeff), zero.gmp_float.get_mpf_t());

        mps_monomial_poly_set_coefficient_f(ctx, mps_poly, i, coeff);
      }
    mpc_clear(coeff);
  }
  ASSERT_EQUAL(mps_monomial_poly_get_precision(ctx, mps_poly), prec);

  // Disable threading.
  // TODO the function is not public:
  // mps_thread_pool_set_concurrency_limit (ctx, NULL, 1);
  // So we have to set thread_safe = false as a workaround.
  // NB: Threads are still created in thread pool when context is created.
  MPS_POLYNOMIAL(mps_poly)->thread_safe = false;
  mps_context_set_input_poly(ctx, MPS_POLYNOMIAL(mps_poly));
  mps_context_set_output_goal(ctx, MPS_OUTPUT_GOAL_APPROXIMATE);
  // MPS_ALGORITHM_STANDARD_MPSOLVE was several times faster
  // in most of our realistic test cases,
  // but sometimes it became incredibly slow (10+ hours instead of 1 minute).
  // Thus, we use MPS_ALGORITHM_SECULAR_GA which seems to be more robust.
  mps_context_select_algorithm(ctx, MPS_ALGORITHM_SECULAR_GA);

  {
    Scoped_Timer mpsolve_timer(timers, "mpsolve");
    mps_mpsolve(ctx);
  }

  mpc_t *mpc_roots = nullptr;
  rdpe_t *radii = nullptr;
  mps_context_get_roots_m(ctx, &mpc_roots, &radii);

  // mps_context_get_degree(ctx) is a degree of zero-deflated polynomial,
  // i.e. number of non-zero roots
  for(int i = 0; i < mps_context_get_degree(ctx); i++)
    {
      El::BigFloat re, im;
      mpf_set(re.gmp_float.get_mpf_t(), mpc_roots[i]->r);
      mpf_set(im.gmp_float.get_mpf_t(), mpc_roots[i]->i);

      // radius = d * 10^l.
      double d;
      long int l;
      rdpe_get_2dl(&d, &l, radii[i]);
      // El::Pow(BigFloat, long) is not implemented, so we use Boost_Float
      El::BigFloat radius = d * to_BigFloat(pow(Boost_Float(10), l));

      result.emplace_back(El::Complex<El::BigFloat>(re, im), radius,
                          mps_context_get_root_status(ctx, i));
    }

  const int num_zero_roots = mps_context_get_zero_roots(ctx);
  if(num_zero_roots > 0)
    {
      const auto status = num_zero_roots == 1 ? MPS_ROOT_STATUS_APPROXIMATED
                                              : MPS_ROOT_STATUS_MULTIPLE;
      result.emplace_back(0, 0, status);
    }

  mpc_vfree(mpc_roots);
  rdpe_vfree(radii);
  mps_monomial_poly_free(ctx, MPS_POLYNOMIAL(mps_poly));
  mps_context_free(ctx);

  if(mps_context_has_errors(ctx))
    RUNTIME_ERROR("MPSolve error: ", mps_context_error_msg(ctx));

  return result;
}

std::vector<El::BigFloat>
find_real_positive_roots_sorted(const Boost_Polynomial &polynomial,
                                Timers &timers)
{
  auto all_roots
    = find_polynomial_roots(to_BigFloat_Vector(polynomial.data()), timers);

  std::vector<El::BigFloat> positive_roots;

  for(const auto &root : all_roots)
    {
      const auto re = root.value.real();
      const auto im = root.value.imag();
      if(re.Sign() <= 0)
        continue;
      // Check that imaginary part is small,
      // using a conservative threshold 2^{-prec/2}
      const auto eps = El::BigFloat(1) >>= El::gmp::Precision() / 2;
      if(Abs(im / re) > eps)
        {
          // TODO check that all complex zeros for conjugate pairs?
          continue;
        }
      ASSERT(MPS_ROOT_STATUS_IS_APPROXIMATED(root.status),
             "Unexpected root status '",
             MPS_ROOT_STATUS_TO_STRING(root.status),
             "' for root=", root.value);

      positive_roots.push_back(re);
    }
  std::sort(positive_roots.begin(), positive_roots.end());

  return positive_roots;
}

std::vector<El::BigFloat>
find_real_positive_minima_sorted(const Boost_Polynomial &polynomial,
                                 Timers &timers)
{
  Scoped_Timer timer(timers, "find_minima");
  std::vector<El::BigFloat> minima;

  // Roots of polynomial derivative
  const auto deriv_roots
    = find_real_positive_roots_sorted(polynomial.prime(), timers);
  if(deriv_roots.empty())
    return minima;

  const size_t num_points = deriv_roots.size();

  // Precompute polynomial values at all extrema
  // We could store only three of them (current, prev, next),
  // but it would complicate the code.
  std::vector<Boost_Float> values;
  {
    values.reserve(num_points);
    for(const auto &x : deriv_roots)
      {
        values.emplace_back(polynomial.evaluate(to_Boost_Float(x)));
      }
  }

  // Value before the first extremum
  const Boost_Float value_zero = polynomial.evaluate(0);
  // Value after the last extremum
  const Boost_Float value_inf
    = polynomial.evaluate(to_Boost_Float(deriv_roots.back() * 2));

  for(size_t i = 0; i < num_points; i++)
    {
      auto &value = values.at(i);
      auto &prev_value = i == 0 ? value_zero : values.at(i - 1);
      auto &next_value = i + 1 == num_points ? value_inf : values.at(i + 1);

      // Check that we have a local minimum indeed
      if(value < prev_value && value < next_value)
        minima.emplace_back(deriv_roots.at(i));
    }

  return minima;
}
