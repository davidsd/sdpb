#include "../factorial.hxx"
#include "pmp_read/pmp_read.hxx"

#include <boost/math/tools/polynomial.hpp>

Boost_Float integral(const Boost_Float &prefactor, const Boost_Float &b,
                     const Boost_Float &x, const int64_t &k);

Boost_Float
rest(const int64_t &m, const Boost_Float &p,
     const std::vector<Boost_Float> &sorted_poles,
     const std::pair<std::vector<Boost_Float>::const_iterator,
                     std::vector<Boost_Float>::const_iterator> &equal_range,
     const int64_t &k);

Boost_Float bilinear_form(
  const Damped_Rational &damped_rational,
  const std::vector<Boost_Float> &sorted_poles,
  const std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                              std::vector<Boost_Float>::const_iterator>>
    &equal_ranges,
  const std::vector<int64_t> &lengths,
  const std::vector<Boost_Float> &products,
  const std::vector<std::vector<Boost_Float>> &integral_matrix,
  const int64_t &m)
{
  Boost_Float result(0);

  size_t index(0);
  for(auto pole(sorted_poles.begin()); pole != sorted_poles.end();)
    {
      const Boost_Float &p(*pole);
      auto &equal_range(equal_ranges.at(index));
      auto &l(lengths.at(index));
      auto &product(products.at(index));

      Boost_Float integral_sum(0);
      auto &integrals(integral_matrix.at(index));
      for(int64_t k = 0; k < l; ++k)
        {
          integral_sum
            += integrals.at(k) * rest(m, p, sorted_poles, equal_range, k);
        }
      result += (pow(p, m) * product) * integral_sum;

      std::advance(pole, l);
      ++index;
    }

  boost::math::tools::polynomial<Boost_Float> numerator(
    pow(boost::math::tools::polynomial<Boost_Float>({0, 1}), m)),
    divisor({1});
  for(auto &pole : sorted_poles)
    {
      divisor *= boost::math::tools::polynomial<Boost_Float>({-pole, 1});
    }
  boost::math::tools::polynomial<Boost_Float> quotient(
    boost::math::tools::quotient_remainder(numerator, divisor).first);

  for(int64_t n = 0; n < int64_t(quotient.size()); ++n)
    {
      result += quotient[n] * factorial(n)
                * pow(-log(damped_rational.base), -1 - n);
    }
  return result * damped_rational.constant;
}
