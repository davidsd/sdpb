// Encoding points for the two-set encoding (pmp2sdp --encodingPoints).
//
// A block's polynomials are defined at its approximation points (pmp2sdp's
// sample points, the nodes of arXiv:2509.14307, as without this option): with
// a reduced prefactor, the block's polynomial program is the degree-d
// polynomial that interpolates R = P / prod(removed poles) there.  The SDP is
// then written at a second set of d + 1 points, the encoding points y_j, with
// scalings w_j: B and c at y_j are w_j times that interpolant at y_j, and the
// bilinear bases are orthonormal for sum_j w_j delta(x - y_j) (even) and
// sum_j w_j y_j delta(x - y_j) (odd).  The polynomial program, and so the
// feasible set and the optimum, are the same as with the default encoding;
// only the SDP's coordinates change (see Dual_Constraint_Group/
// two_set_encoding.cxx for the second step).
//
// Rules for the encoding points of a block with m kept (reduced-prefactor)
// poles p_i <= 0 and d + 1 approximation points, with the scale
// a = A * median_i |p_i| (a = 10 A without poles; or a per-block a):
//
//   functionalMatched: the quadrature rule for the measure
//     dx / (sqrt(x) g(x)) on [0, infinity) whose points are the "midpoints"
//     of the counting function N_g.  With
//       g(x)   = prod_i (x + |p_i|) (x + a)^e,   e = d + 1 - m,
//       N_g(x) = (2/pi) [ sum_i atan sqrt(x/|p_i|) + e atan sqrt(x/a) ],
//       P_g(x) = (g(x)/pi) [ sum_i sqrt|p_i| / (x + |p_i|)
//                            + e sqrt(a) / (x + a) ],
//     the points solve N_g(y_k) = k - 1/2 (k = 1..d+1) and the scalings are
//     w_k = 1 / P_g(y_k).  Then sum_k w_k p(y_k) = int_0^inf p dx /
//     (sqrt(x) g(x)) for every polynomial p of degree <= d.  A pole at x = 0
//     has no finite home measure; such poles are counted in (x + a)^e.
//   chebyshev: the Chebyshev points of t = x / (x + a),
//     y_k = a tan^2((2k - 1) pi / (4 (d + 1))), with the scalings
//     w_k = 1 / (prod_i (y_k - p_i) (y_k + a)^(d - m)) = 1 / P_ref(y_k).
//
// Both put most points far beyond the approximation points, on the uncut
// half-line, so the cone is exact.

#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/math/constants/constants.hpp>

#include <algorithm>
#include <cmath>

namespace
{
  std::vector<Boost_Float> pole_magnitudes(const Damped_Rational &prefactor)
  {
    std::vector<Boost_Float> magnitudes;
    for(const auto &p : prefactor.poles)
      {
        ASSERT(p <= 0, "All poles must be less or equal to zero!",
               DEBUG_STRING(p), DEBUG_STRING(prefactor));
        magnitudes.push_back(-p);
      }
    std::sort(magnitudes.begin(), magnitudes.end());
    return magnitudes;
  }
}

El::BigFloat encoding_points_scale(const Damped_Rational &reduced_prefactor,
                                   const El::BigFloat &scale_factor)
{
  ASSERT(scale_factor > 0, DEBUG_STRING(scale_factor));
  const auto magnitudes = pole_magnitudes(reduced_prefactor);
  const Boost_Float a
    = to_Boost_Float(scale_factor)
      * (magnitudes.empty() ? Boost_Float(10)
                            : magnitudes.at(magnitudes.size() / 2));
  ASSERT(a > 0, "the encoding points' scale a = factor * median |pole| must "
                "be positive (all kept poles at x = 0?); give it per block "
                "with --encodingScaleFile",
         DEBUG_STRING(reduced_prefactor));
  return to_BigFloat(a);
}

void functional_matched_encoding_points(
  const Damped_Rational &reduced_prefactor, const size_t num_points,
  const El::BigFloat &scale, std::vector<El::BigFloat> &points,
  std::vector<El::BigFloat> &scalings)
{
  ASSERT(num_points > 0);
  ASSERT(scale > 0, DEBUG_STRING(scale));
  const Boost_Float a = to_Boost_Float(scale);
  const Boost_Float pi = boost::math::constants::pi<Boost_Float>();

  std::vector<Boost_Float> q;
  int64_t e = static_cast<int64_t>(num_points);
  for(const auto &m : pole_magnitudes(reduced_prefactor))
    if(m > 0)
      {
        q.push_back(m);
        --e;
      }
  ASSERT(e >= 0,
         "--encodingPoints functionalMatched: more nonzero kept poles than "
         "points",
         DEBUG_STRING(q.size()), DEBUG_STRING(num_points));

  // N_g and its derivative, in working precision
  const auto N = [&](const Boost_Float &x) {
    Boost_Float s = e * atan(sqrt(x / a));
    for(const auto &qi : q)
      s += atan(sqrt(x / qi));
    return Boost_Float(2 * s / pi);
  };
  const auto dN = [&](const Boost_Float &x) {
    Boost_Float s = e * sqrt(a) / (x + a);
    for(const auto &qi : q)
      s += sqrt(qi) / (x + qi);
    return Boost_Float(s / (pi * sqrt(x)));
  };
  const auto P_g = [&](const Boost_Float &x) {
    Boost_Float g = pow(x + a, e);
    Boost_Float s = e * sqrt(a) / (x + a);
    for(const auto &qi : q)
      {
        g *= x + qi;
        s += sqrt(qi) / (x + qi);
      }
    return Boost_Float(g * s / pi);
  };
  // A start in double precision, by bisection in log(x)
  std::vector<double> qd;
  for(const auto &qi : q)
    qd.push_back(static_cast<double>(log(qi)));
  const double log_a = static_cast<double>(log(a));
  const auto N_double = [&](const double log_x) {
    double s = e * std::atan(std::exp(0.5 * (log_x - log_a)));
    for(const auto &lq : qd)
      s += std::atan(std::exp(0.5 * (log_x - lq)));
    return 2 * s / M_PI;
  };

  const auto digits2 = Boost_Float::default_precision() * std::log2(10.0);
  points.resize(num_points);
  scalings.resize(num_points);
  for(size_t k = 0; k < num_points; ++k)
    {
      const double target = k + 0.5;
      double lo = -1400, hi = 1400;
      for(int it = 0; it < 200; ++it)
        {
          const double mid = 0.5 * (lo + hi);
          (N_double(mid) < target ? lo : hi) = mid;
        }
      Boost_Float x = exp(Boost_Float(0.5 * (lo + hi)));
      const Boost_Float target_bf = Boost_Float(2 * k + 1) / 2;
      bool converged = false;
      for(int it = 0; it < 64 && !converged; ++it)
        {
          const Boost_Float dx = (N(x) - target_bf) / dN(x);
          x -= dx;
          ASSERT(x > 0, "--encodingPoints functionalMatched: Newton left "
                        "the half-line",
                 DEBUG_STRING(k));
          // Converged when the step is at the working precision, or when
          // N(x) is (the far points' relative accuracy is limited by
          // x N'(x), which is small there).
          converged = abs(dx) <= x * pow(Boost_Float(2), -digits2 + 8)
                      || abs(N(x) - target_bf)
                           <= num_points * pow(Boost_Float(2), -digits2 + 8);
        }
      ASSERT(converged,
             "--encodingPoints functionalMatched: Newton did not converge",
             DEBUG_STRING(k));
      points.at(k) = to_BigFloat(x);
      scalings.at(k) = to_BigFloat(Boost_Float(1 / P_g(x)));
    }
  for(size_t k = 1; k < num_points; ++k)
    ASSERT(points.at(k) > points.at(k - 1), DEBUG_STRING(k));
}

void chebyshev_encoding_points(const Damped_Rational &reduced_prefactor,
                               const size_t num_points,
                               const El::BigFloat &scale,
                               std::vector<El::BigFloat> &points,
                               std::vector<El::BigFloat> &scalings)
{
  ASSERT(num_points > 0);
  ASSERT(scale > 0, DEBUG_STRING(scale));
  const Boost_Float a = to_Boost_Float(scale);
  const Boost_Float pi = boost::math::constants::pi<Boost_Float>();
  const int64_t degree = static_cast<int64_t>(num_points) - 1;
  const int64_t exponent
    = degree - static_cast<int64_t>(reduced_prefactor.poles.size());
  points.resize(num_points);
  scalings.resize(num_points);
  for(size_t i = 0; i < num_points; ++i)
    {
      const Boost_Float u = (2 * i + 1) * pi / (4 * num_points);
      const Boost_Float t = tan(u);
      const Boost_Float x = a * t * t;
      Boost_Float denominator = pow(x + a, exponent);
      for(const auto &p : reduced_prefactor.poles)
        denominator *= x - p;
      points.at(i) = to_BigFloat(x);
      scalings.at(i) = to_BigFloat(Boost_Float(1 / denominator));
    }
}

void set_encoding_points(Polynomial_Vector_Matrix &m,
                         std::vector<El::BigFloat> points,
                         std::vector<El::BigFloat> scalings)
{
  const size_t num_points = m.sample_points.size();
  ASSERT_EQUAL(points.size(), num_points,
               "the number of encoding points must be the number of "
               "approximation points");
  ASSERT_EQUAL(scalings.size(), num_points);
  for(size_t k = 0; k < num_points; ++k)
    {
      ASSERT(points.at(k) >= 0, "encoding points must be on x >= 0",
             DEBUG_STRING(k), DEBUG_STRING(points.at(k)));
      ASSERT(scalings.at(k) > 0, "encoding scalings must be positive",
             DEBUG_STRING(k), DEBUG_STRING(scalings.at(k)));
    }
  {
    auto sorted = points;
    std::sort(sorted.begin(), sorted.end());
    for(size_t k = 1; k < num_points; ++k)
      ASSERT(sorted.at(k) != sorted.at(k - 1),
             "encoding points must be distinct", DEBUG_STRING(k));
  }
  m.encoding_points = std::move(points);
  m.encoding_scalings = std::move(scalings);
}
