#include "../Damped_Rational.hxx"

#include <El.hpp>

#include <boost/math/constants/constants.hpp>

#include <vector>

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational)
{
  std::vector<Boost_Float> result;
  result.reserve(points.size());

  // Evaluate damped_rational at given points.
  // If a point x is too close to a pole p,
  // we replace (x - p) with min_pole_distance.

  // TODO choose better threshold and/or make it configurable
  const Boost_Float min_pole_distance("1e-16");

  for(auto &point : points)
    {
      Boost_Float numerator(damped_rational.constant
                            * pow(damped_rational.base, point));
      Boost_Float denominator(1);
      for(auto &pole : damped_rational.poles)
        {
          Boost_Float delta = point - pole;

          // Regularize value at small denominator
          if(abs(delta) < min_pole_distance)
            {
              delta = min_pole_distance;
              if(sign(delta) < 0)
                delta = -delta;
            }

          denominator *= delta;
        }
      result.emplace_back(numerator / denominator);
    }
  return result;
}
