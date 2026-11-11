// Interpolation conditions at x = infinity for a reduced prefactor
// (pmp2sdp --reducedPrefactorNodesAtInfinity M).
//
// With a reduced prefactor chi_red (fewer poles than the prefactor chi),
// the rows of a block at the sample points x_0 < ... < x_d are
//
//   V_k = chi(x_k) P(x_k) = chi_red(x_k) f(x_k),
//   f(x) = (chi/chi_red)(x) P(x) = (c/c_red) P(x) Pi_red(x) / Pi(x),
//
// where Pi(x) = prod_{poles p} (x - p), so the SDP sees the degree-d
// polynomial Q that interpolates f at all d+1 sample points. Q is not exact
// at x = infinity: its leading coefficient differs from that of f.
//
// With M conditions at infinity we use instead the degree-d polynomial Qinf
// that interpolates f at the d+1-M smallest sample points and whose top M
// coefficients (of x^d, ..., x^{d-M+1}) equal those of f's expansion at
// infinity, i.e. of the polynomial part of (c/c_red) P Pi_red / Pi. The
// sample points, bilinear basis and the rows at the d+1-M smallest points are
// unchanged; the rows at the M largest points become chi_red(x_i) Qinf(x_i).
//
// Construction: with F the d+1-M smallest points, L_F[f] their Lagrange
// interpolant (degree <= d-M) and omega_F(x) = prod_{k in F}(x - x_k),
//
//   Qinf = L_F[f] + omega_F R,
//
// where R (degree M-1) is fixed by the top M coefficients, which L_F[f] does
// not touch. The arithmetic is done in MPFR at twice the working precision.

#include "pmp2sdp/Dual_Constraint_Group.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Boost_Float.hxx"

#include <mpfr.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace
{
  // Fixed-size array of MPFR numbers, all with the same precision.
  class Mpfr_Array
  {
  public:
    Mpfr_Array(const size_t size, const mpfr_prec_t precision) : data(size)
    {
      for(auto &x : data)
        {
          mpfr_init2(&x, precision);
          mpfr_set_zero(&x, 1);
        }
    }
    ~Mpfr_Array()
    {
      for(auto &x : data)
        mpfr_clear(&x);
    }
    Mpfr_Array(const Mpfr_Array &) = delete;
    Mpfr_Array &operator=(const Mpfr_Array &) = delete;

    mpfr_ptr operator[](const size_t index) { return &data.at(index); }
    [[nodiscard]] size_t size() const { return data.size(); }

  private:
    std::vector<__mpfr_struct> data;
  };

  void set(mpfr_ptr result, const Boost_Float &value)
  {
    mpfr_set(result, value.backend().data(), MPFR_RNDN);
  }
  void set(mpfr_ptr result, const El::BigFloat &value)
  {
    set(result, to_Boost_Float(value));
  }
  El::BigFloat mpfr_to_BigFloat(mpfr_srcptr value)
  {
    Boost_Float result;
    mpfr_set_prec(result.backend().data(), El::gmp::Precision());
    mpfr_set(result.backend().data(), value, MPFR_RNDN);
    return ::to_BigFloat(result);
  }

  double log10_abs(mpfr_srcptr value)
  {
    long exponent;
    const double mantissa = mpfr_get_d_2exp(&exponent, value, MPFR_RNDN);
    return std::log10(std::abs(mantissa))
           + static_cast<double>(exponent) * std::log10(2.0);
  }

  // Top coefficients of prod_k (x - roots[k]) (monic):
  // result[j] = coefficient of x^{n-j}, j = 0..result.size()-1.
  void top_coefficients_from_roots(Mpfr_Array &roots, Mpfr_Array &result,
                                   const mpfr_prec_t precision)
  {
    Mpfr_Array tmp(1, precision);
    for(size_t j = 0; j < result.size(); ++j)
      mpfr_set_ui(result[j], j == 0 ? 1 : 0, MPFR_RNDN);
    for(size_t k = 0; k < roots.size(); ++k)
      // multiply by (x - root): c_j <- c_j - root * c_{j-1}
      for(size_t j = result.size() - 1; j >= 1; --j)
        {
          mpfr_mul(tmp[0], roots[k], result[j - 1], MPFR_RNDN);
          mpfr_sub(result[j], result[j], tmp[0], MPFR_RNDN);
        }
  }
}

void apply_nodes_at_infinity(Dual_Constraint_Group &group,
                             const Polynomial_Vector_Matrix &m,
                             const size_t num_nodes_at_infinity)
{
  if(num_nodes_at_infinity == 0)
    return;

  // Same poles: chi/chi_red is a constant, f is a polynomial of degree <= d,
  // and the interpolant is already exact.
  {
    auto poles = m.prefactor.poles;
    auto reduced_poles = m.reduced_prefactor.poles;
    std::sort(poles.begin(), poles.end());
    std::sort(reduced_poles.begin(), reduced_poles.end());
    if(poles == reduced_poles)
      return;
  }
  ASSERT(m.prefactor.base == m.reduced_prefactor.base,
         "--reducedPrefactorNodesAtInfinity requires prefactor and "
         "reducedPrefactor with the same base",
         DEBUG_STRING(m.prefactor), DEBUG_STRING(m.reduced_prefactor));

  const auto &polys = m.polynomials;
  const size_t num_points = group.num_points;
  const size_t dim = group.dim;
  const size_t vector_dim = polys(0, 0).size();
  const size_t M = std::min(num_nodes_at_infinity, num_points);
  const size_t num_finite = num_points - M;
  const size_t num_poles = m.prefactor.poles.size();
  const size_t num_reduced_poles = m.reduced_prefactor.poles.size();

  int64_t max_degree = 0;
  for(size_t r = 0; r < polys.Height(); ++r)
    for(size_t c = 0; c < polys.Width(); ++c)
      for(const auto &poly : polys(r, c))
        max_degree = std::max(max_degree, poly.degree());
  // d = num_points - 1 = max_degree - (num_poles - num_reduced_poles)
  ASSERT_EQUAL(static_cast<int64_t>(num_points),
               max_degree + 1 + static_cast<int64_t>(num_reduced_poles)
                 - static_cast<int64_t>(num_poles));
  ASSERT_EQUAL(m.sample_points.size(), num_points);
  ASSERT_EQUAL(m.reduced_sample_scalings.size(), num_points);

  const mpfr_prec_t precision = 2 * El::gmp::Precision() + 64;

  // Sample points in increasing order: F = the num_finite smallest,
  // I = the M largest (moved to infinity).
  std::vector<size_t> order(num_points);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&m](size_t a, size_t b) {
    return m.sample_points.at(a) < m.sample_points.at(b);
  });

  Mpfr_Array x(num_points, precision);
  Mpfr_Array scaling(num_points, precision);
  for(size_t k = 0; k < num_points; ++k)
    {
      set(x[k], m.sample_points.at(k));
      set(scaling[k], m.reduced_sample_scalings.at(k));
    }

  Mpfr_Array tmp(2, precision);

  // weights[i * num_finite + a] = l_a(x_I[i]), the Lagrange basis polynomial
  // of the a-th point of F, evaluated at the i-th point of I.
  // omega[i] = omega_F(x_I[i]).
  Mpfr_Array weights(M * num_finite, precision);
  Mpfr_Array omega(M, precision);
  for(size_t i = 0; i < M; ++i)
    {
      const auto xi = x[order.at(num_finite + i)];
      mpfr_set_ui(omega[i], 1, MPFR_RNDN);
      for(size_t a = 0; a < num_finite; ++a)
        {
          mpfr_sub(tmp[0], xi, x[order.at(a)], MPFR_RNDN);
          mpfr_mul(omega[i], omega[i], tmp[0], MPFR_RNDN);
        }
      for(size_t a = 0; a < num_finite; ++a)
        {
          const auto xa = x[order.at(a)];
          auto w = weights[i * num_finite + a];
          mpfr_set_ui(w, 1, MPFR_RNDN);
          for(size_t b = 0; b < num_finite; ++b)
            {
              if(b == a)
                continue;
              const auto xb = x[order.at(b)];
              mpfr_sub(tmp[0], xi, xb, MPFR_RNDN);
              mpfr_sub(tmp[1], xa, xb, MPFR_RNDN);
              ASSERT(!mpfr_zero_p(tmp[1]), "Repeated sample point");
              mpfr_div(tmp[0], tmp[0], tmp[1], MPFR_RNDN);
              mpfr_mul(w, w, tmp[0], MPFR_RNDN);
            }
        }
    }

  // Top M coefficients of omega_F, Pi and Pi_red (all monic).
  Mpfr_Array omega_top(M, precision), pi_top(M, precision),
    pi_reduced_top(M, precision);
  {
    Mpfr_Array roots(num_finite, precision);
    for(size_t a = 0; a < num_finite; ++a)
      mpfr_set(roots[a], x[order.at(a)], MPFR_RNDN);
    top_coefficients_from_roots(roots, omega_top, precision);
  }
  {
    Mpfr_Array roots(num_poles, precision);
    for(size_t k = 0; k < num_poles; ++k)
      set(roots[k], m.prefactor.poles.at(k));
    top_coefficients_from_roots(roots, pi_top, precision);
  }
  {
    Mpfr_Array roots(num_reduced_poles, precision);
    for(size_t k = 0; k < num_reduced_poles; ++k)
      set(roots[k], m.reduced_prefactor.poles.at(k));
    top_coefficients_from_roots(roots, pi_reduced_top, precision);
  }
  // c / c_red
  Mpfr_Array ratio(1, precision);
  {
    Mpfr_Array c_red(1, precision);
    set(ratio[0], m.prefactor.constant);
    set(c_red[0], m.reduced_prefactor.constant);
    mpfr_div(ratio[0], ratio[0], c_red[0], MPFR_RNDN);
  }

  // The moved rows are computed from f(x_k) = row / reducedSampleScaling_k
  // and from q, which comes from the prefactors' poles and constants. The two
  // agree only if sampleScaling_k / reducedSampleScaling_k
  // = (c/c_red) Pi_red(x_k) / Pi(x_k), as for pmp2sdp's default scalings.
  // Explicit scalings that differ are refused.
  ASSERT_EQUAL(m.sample_scalings.size(), num_points);
  {
    Mpfr_Array expected(1, precision), actual(1, precision);
    const double max_log2_error = -0.5 * El::gmp::Precision();
    for(size_t k = 0; k < num_points; ++k)
      {
        mpfr_set(expected[0], ratio[0], MPFR_RNDN);
        for(size_t j = 0; j < num_reduced_poles; ++j)
          {
            set(tmp[0], m.reduced_prefactor.poles.at(j));
            mpfr_sub(tmp[0], x[k], tmp[0], MPFR_RNDN);
            mpfr_mul(expected[0], expected[0], tmp[0], MPFR_RNDN);
          }
        bool at_pole = false;
        for(size_t j = 0; j < num_poles; ++j)
          {
            set(tmp[0], m.prefactor.poles.at(j));
            mpfr_sub(tmp[0], x[k], tmp[0], MPFR_RNDN);
            at_pole = at_pole || mpfr_zero_p(tmp[0]);
            mpfr_div(expected[0], expected[0], tmp[0], MPFR_RNDN);
          }
        set(actual[0], m.sample_scalings.at(k));
        set(tmp[1], m.reduced_sample_scalings.at(k));
        mpfr_div(actual[0], actual[0], tmp[1], MPFR_RNDN);
        mpfr_sub(tmp[0], actual[0], expected[0], MPFR_RNDN);
        const bool agree
          = !at_pole
            && (mpfr_zero_p(tmp[0])
                || (!mpfr_zero_p(expected[0])
                    && mpfr_get_exp(tmp[0]) - mpfr_get_exp(expected[0])
                         < max_log2_error));
        ASSERT(agree,
               "--reducedPrefactorNodesAtInfinity requires sampleScalings / "
               "reducedSampleScalings = the ratio of the prefactors at each "
               "sample point, as pmp2sdp's default scalings are; it differs "
               "at sample point ",
               k, DEBUG_STRING(m.prefactor), DEBUG_STRING(m.reduced_prefactor));
      }
  }

  Mpfr_Array P_top(M, precision), N_top(M, precision), q(M, precision),
    R(M, precision), f(num_finite, precision), value(1, precision),
    magnitude(1, precision);
  double worst_cancellation = 0;

  size_t p0 = 0; // first row of the current (r, c) entry
  for(size_t c = 0; c < dim; ++c)
    for(size_t r = 0; r <= c; ++r, p0 += num_points)
      for(size_t n = 0; n < vector_dim; ++n)
        {
          // Rows are  c_p = V and B_{p, n-1} = -V
          const auto row = [&group, n](const size_t p) -> El::BigFloat {
            return n == 0 ? group.constraint_constants.at(p)
                          : -group.constraint_matrix.Get(p, n - 1);
          };
          const auto set_row = [&group, n](const size_t p,
                                           const El::BigFloat &v) {
            if(n == 0)
              group.constraint_constants.at(p) = v;
            else
              group.constraint_matrix.Set(p, n - 1, -v);
          };

          // q[j] = coefficient of x^{d-j} in the polynomial part of
          // (c/c_red) P Pi_red / Pi
          const auto &coefficients = polys(r, c).at(n).coefficients;
          for(size_t l = 0; l < M; ++l)
            {
              const int64_t power = max_degree - static_cast<int64_t>(l);
              if(power >= 0
                 && static_cast<size_t>(power) < coefficients.size())
                set(P_top[l], coefficients.at(power));
              else
                mpfr_set_zero(P_top[l], 1);
            }
          for(size_t j = 0; j < M; ++j)
            {
              mpfr_set_zero(N_top[j], 1);
              for(size_t l = 0; l <= j; ++l)
                {
                  mpfr_mul(tmp[0], P_top[l], pi_reduced_top[j - l],
                           MPFR_RNDN);
                  mpfr_add(N_top[j], N_top[j], tmp[0], MPFR_RNDN);
                }
              mpfr_set(q[j], N_top[j], MPFR_RNDN);
              for(size_t l = 1; l <= j; ++l)
                {
                  mpfr_mul(tmp[0], q[j - l], pi_top[l], MPFR_RNDN);
                  mpfr_sub(q[j], q[j], tmp[0], MPFR_RNDN);
                }
            }
          for(size_t j = 0; j < M; ++j)
            mpfr_mul(q[j], q[j], ratio[0], MPFR_RNDN);

          // R[j] = coefficient of x^{M-1-j} in R, from the top M
          // coefficients of omega_F R = q
          for(size_t j = 0; j < M; ++j)
            {
              mpfr_set(R[j], q[j], MPFR_RNDN);
              for(size_t l = 1; l <= j; ++l)
                {
                  mpfr_mul(tmp[0], R[j - l], omega_top[l], MPFR_RNDN);
                  mpfr_sub(R[j], R[j], tmp[0], MPFR_RNDN);
                }
            }

          // f at the points of F, from the SDP's own rows
          for(size_t a = 0; a < num_finite; ++a)
            {
              const size_t k = order.at(a);
              set(f[a], row(p0 + k));
              mpfr_div(f[a], f[a], scaling[k], MPFR_RNDN);
            }

          for(size_t i = 0; i < M; ++i)
            {
              const size_t k = order.at(num_finite + i);
              const auto xi = x[k];
              mpfr_set_zero(value[0], 1);
              mpfr_set_zero(magnitude[0], 1);
              for(size_t a = 0; a < num_finite; ++a)
                {
                  mpfr_mul(tmp[0], weights[i * num_finite + a], f[a],
                           MPFR_RNDN);
                  mpfr_add(value[0], value[0], tmp[0], MPFR_RNDN);
                  mpfr_abs(tmp[0], tmp[0], MPFR_RNDN);
                  mpfr_add(magnitude[0], magnitude[0], tmp[0], MPFR_RNDN);
                }
              // R(x_i) by Horner
              mpfr_set(tmp[1], R[0], MPFR_RNDN);
              for(size_t j = 1; j < M; ++j)
                {
                  mpfr_mul(tmp[1], tmp[1], xi, MPFR_RNDN);
                  mpfr_add(tmp[1], tmp[1], R[j], MPFR_RNDN);
                }
              mpfr_mul(tmp[1], tmp[1], omega[i], MPFR_RNDN);
              mpfr_add(value[0], value[0], tmp[1], MPFR_RNDN);
              mpfr_abs(tmp[1], tmp[1], MPFR_RNDN);
              mpfr_add(magnitude[0], magnitude[0], tmp[1], MPFR_RNDN);

              // Digits lost to cancellation in Qinf(x_i)
              if(!mpfr_zero_p(magnitude[0]))
                {
                  const double lost
                    = mpfr_zero_p(value[0])
                        ? static_cast<double>(precision) * std::log10(2.0)
                        : log10_abs(magnitude[0]) - log10_abs(value[0]);
                  worst_cancellation = std::max(worst_cancellation, lost);
                }

              mpfr_mul(value[0], value[0], scaling[k], MPFR_RNDN);
              set_row(p0 + k, mpfr_to_BigFloat(value[0]));
            }
        }
  group.nodes_at_infinity_digits_lost = worst_cancellation;
}
