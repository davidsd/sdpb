#include "../Damped_Rational.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"

#include <El.hpp>

#include <boost/math/constants/constants.hpp>

#include <vector>

// TODO add comments describing sampling algorithm and its rationale.

// TODO choose better threshold and/or make it configurable
#define SMALL_POLE_THRESHOLD Boost_Float("1e-10")

namespace
{
  // Truncate x to [-1,1] range (to compensate rounding errors) and evaluate acos(x)
  Boost_Float acos_safe(const Boost_Float &x)
  {
    if(abs(x) > 1)
      {
        // Print warning if x is too far from [-1,1].
        // If delta is ~ rounding error (i.e. 1/2^{binary precision}), it's OK.
        // If delta is much larger, it's suspicious and we print a warning.
        const Boost_Float eps
          = to_Boost_Float(El::BigFloat(1) >>= El::gmp::Precision() / 2);
        if(abs(x) > 1 + eps)
          {
            PRINT_WARNING("acos argument lies outside of [-1,1] range and "
                          "will be truncated",
                          DEBUG_STRING(x), DEBUG_STRING(abs(x) - 1));
          }

        // Set x to -1 or 1
        return acos(Boost_Float(sign(x)));
      }
    return acos(x);
  }

  // b equation and its first derivative
  std::tuple<Boost_Float, Boost_Float>
  b_equation(const size_t &num_points, const Damped_Rational &prefactor,
             const Boost_Float &b)
  {
    // bEquation = Sum[1-Sqrt[-p/(bVar-p)],{p,poles}]-1/2 bVar Log[base]-numSamplePoints
    // bEquation' = -(Log[base]/2) + Sum[(p/(-bVar + p))^(3/2)/(2 (-p)), {p, poles}]
    Boost_Float eq = 0;
    Boost_Float eq_deriv = 0;
    for(const auto &p : prefactor.poles)
      {
        eq += 1 - sqrt(-p / (b - p));
        eq_deriv += sqrt(-p) / pow(sqrt(b - p), 3) / 2;
      }
    eq += -b * log(prefactor.base) / 2 - num_points;
    eq_deriv += -log(prefactor.base) / 2;

    ASSERT(!isnan(eq) && !isnan(eq_deriv) && !isinf(eq) && !isinf(eq_deriv),
           DEBUG_STRING(num_points), DEBUG_STRING(prefactor), DEBUG_STRING(b),
           DEBUG_STRING(eq), DEBUG_STRING(eq_deriv));
    return {eq, eq_deriv};
  }

  Boost_Float
  find_b(const size_t &num_points, const Damped_Rational &prefactor)
  {
    const auto F = [&](const Boost_Float &b) {
      return b_equation(num_points, prefactor, b);
    };

    // step away from zero to avoid 0/0 in b_equation at p=0
    const Boost_Float min = SMALL_POLE_THRESHOLD;
    const Boost_Float max = -(2 * num_points / log(prefactor.base));
    ASSERT(min <= max, DEBUG_STRING(min), DEBUG_STRING(max),
           DEBUG_STRING(num_points), DEBUG_STRING(prefactor));

    const Boost_Float guess = (min + max) / 2;
    // TODO: can we use full precision?
    const auto digits2 = El::gmp::Precision() / 2;
    return boost::math::tools::newton_raphson_iterate(F, guess, min, max,
                                                      digits2);
  }

  // integratedDensity (as function of z) and its derivative
  // integratedDensity = Sum[
  //    1 /\[Pi] ( ArcCos[1-(2z(b-p))/(b (z-p))]
  //    - Sqrt[-p/(b-p)] ArcCos[1-(2 z)/b])
  //    ,{p,poles}
  //    ]
  //    -Log[base]/\[Pi] (Sqrt[(b-z) z]+ b/2 ArcCos[1-(2 z)/b]);
  //
  // integratedDensity'[z] == Sqrt[b - z]/(\[Pi] Sqrt[z]) (Sum[Sqrt[-p]/(Sqrt[b - p] (z - p)), {p, poles}] - Log[base])
  std::tuple<Boost_Float, Boost_Float>
  integrated_density(const Damped_Rational &prefactor, const Boost_Float &b,
                     const Boost_Float &z)
  {
    ASSERT(z <= b, DEBUG_STRING(z), DEBUG_STRING(b), DEBUG_STRING(prefactor));
    const auto pi = boost::math::constants::pi<Boost_Float>();

    // NB: we cannot initialize one_div_pi as a global constant,
    // because we don't know required MPFR precision at that moment.
    // Thus, we do it inside the function.
    const Boost_Float one_div_pi =
#if(BOOST_VERSION >= 107200)
      boost::math::constants::one_div_pi<Boost_Float>();
#else
      Boost_Float(1) / pi;
#endif

    Boost_Float density = 0;
    Boost_Float density_deriv = 0;
    for(const auto &p : prefactor.poles)
      {
        // TODO there is a square-root singularity at p = 0,
        // which leads to precision loss. But it is probably not important for us.
        density += one_div_pi
                   * (acos_safe(1 - (2 * z * (b - p)) / (b * (z - p)))
                      - sqrt(-p / (b - p)) * acos_safe(1 - (2 * z) / b));
        density_deriv
          += sqrt(-p) / (sqrt(b - p) * (z - p)) * sqrt(b - z) / (pi * sqrt(z));
      }
    density += -log(prefactor.base) / pi
               * (sqrt((b - z) * z) + b / 2 * acos_safe(1 - (2 * z) / b));
    density_deriv += -log(prefactor.base) * sqrt(b - z) / (pi * sqrt(z));

    ASSERT(!isnan(density) && !isnan(density_deriv) && !isinf(density)
             && !isinf(density_deriv),
           DEBUG_STRING(prefactor), DEBUG_STRING(b), DEBUG_STRING(z),
           DEBUG_STRING(density), DEBUG_STRING(density_deriv));
    return {density, density_deriv};
  }

  // Solve equation: integrated_density(z) == n + 1/2
  // and write root to output[n]
  // for n in [n_min..num_points)
  void find_bohr_sommerfeld_roots(const Damped_Rational &prefactor,
                                  const size_t &n_min,
                                  std::vector<Boost_Float> &output)
  {
    const auto num_points = output.size();
    ASSERT(num_points > 0);
    ASSERT(n_min <= num_points, DEBUG_STRING(n_min), DEBUG_STRING(num_points));

    // No more points needed
    if(n_min == num_points)
      return;

    const auto b = find_b(num_points, prefactor);
    ASSERT(b > 0, DEBUG_STRING(b), DEBUG_STRING(num_points),
           DEBUG_STRING(prefactor));

    Boost_Float min = SMALL_POLE_THRESHOLD;
    const Boost_Float max = b;
    const auto digits2 = El::gmp::Precision() / 2;

    for(size_t n = n_min; n < num_points; ++n)
      {
        ASSERT(min <= max);
        // This guess works well for evenly spaced points
        Boost_Float guess = min + (max - min) / (num_points - n + 1);
        // extra safety measure, for the (unlikely) case of rounding errors:
        guess = std::max(guess, min);
        guess = std::min(guess, max);

        auto F
          = [&](const Boost_Float &z) -> std::tuple<Boost_Float, Boost_Float> {
          const auto [f, f_deriv] = integrated_density(prefactor, b, z);
          return {f - n - 0.5, f_deriv};
        };

        output.at(n) = boost::math::tools::newton_raphson_iterate(
          F, guess, min, max, digits2);
        min = output.at(n);
      }
  }
}

std::vector<Boost_Float>
sample_points(const size_t &num_points, const Damped_Rational &prefactor)
{
  if(num_points == 1)
    {
      if(!prefactor.poles.empty())
        {
          PRINT_WARNING("Prefactor for a constant constraint has poles: ",
                        DEBUG_STRING(prefactor));
        }
      // For constant constraints (degree-0 polynomials),
      // we can choose any sample point, so we choose 0.
      // Usually there will be a constant prefactor = 1,
      // or a default prefactor exp(-x),
      // both leading to trivial sample_scalings={1}.
      return {0};
    }

  // Calculate analytic sample points

  std::vector<Boost_Float> points(num_points);

  // For base>=1, integrals diverge and we cannot get meaningful answer.
  ASSERT(prefactor.base > 0 && prefactor.base < 1, DEBUG_STRING(prefactor));

  // Number of small points
  size_t num_small_points = std::count_if(
    prefactor.poles.begin(), prefactor.poles.end(), [&](const auto &p) {
      ASSERT(p <= 0, "All poles must be less or equal to zero!",
             DEBUG_STRING(p));
      return boost::multiprecision::abs(p) <= SMALL_POLE_THRESHOLD;
    });
  num_small_points = std::min(num_small_points, num_points);

  // Fill points[num_small_points..num_points)
  find_bohr_sommerfeld_roots(prefactor, num_small_points, points);

  // Fill points[0..num_small_points) with evenly spaced values
  // from 0 to the first Bohr-Sommerfeld root (or to b)
  auto small_point_end = num_small_points == num_points
                           ? find_b(num_points, prefactor)
                           : points.at(num_small_points);

  ASSERT(small_point_end > 0, "Cannot sample points near zero",
         DEBUG_STRING(num_points), DEBUG_STRING(prefactor),
         DEBUG_STRING(num_small_points), DEBUG_STRING(small_point_end));
  for(size_t i = 0; i < num_small_points; ++i)
    {
      points.at(i) = small_point_end * i / num_small_points;
    }

  // Sanity check
  for(size_t i = 1; i < num_points; ++i)
    {
      ASSERT(points.at(i) > points.at(i - 1), DEBUG_STRING(i),
             DEBUG_STRING(points.at(i)), DEBUG_STRING(points.at(i - 1)),
             DEBUG_STRING(num_points), DEBUG_STRING(num_small_points),
             DEBUG_STRING(prefactor));
    }

  return points;
}
