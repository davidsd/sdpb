#include "../Damped_Rational.hxx"

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
      result.emplace_back(damped_rational.evaluate(point, min_pole_distance));
    }
  return result;
}
