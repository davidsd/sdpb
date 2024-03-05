#include "../Damped_Rational.hxx"

#include <El.hpp>

#include <boost/math/constants/constants.hpp>

#include <vector>

// Rescaled Laguerre
std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational)
{
  std::vector<Boost_Float> result;
  result.reserve(points.size());
  for(auto &point : points)
    {
      Boost_Float numerator(damped_rational.constant
                            * pow(damped_rational.base, point));
      Boost_Float denominator(1);
      for(auto &pole : damped_rational.poles)
        {
          denominator *= (point - pole);
        }
      result.emplace_back(numerator / denominator);
    }
  return result;
}
