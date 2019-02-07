#include "factorial.hxx"

#include <boost/math/special_functions/expint.hpp>

// integral(0) = b^x tgamma(0,x*log(b)) = -(b^x) expint(0, -x log(b))
// integral(k) = d(integral(k-1),x)/k
//
// In closed form
//
// integral(k) = -(b^x expint(0, -x log(b)) * (log(b))^k
//                + Sum(i={1..=k}) {(-1)^k * (i-1)! (-log(b))^(k-i) / x^i})/k!

Boost_Float
integral(const Boost_Float &b, const Boost_Float &x, const int64_t &k)
{
  Boost_Float result(-boost::math::expint(-x * log(b)) * pow(b, x)
                     * pow(log(b), k));
  for(int64_t i = 1; i <= k; ++i)
    {
      result += (k % 2 == 0 ? 1 : -1) * factorial(i - 1) * pow(-log(b), k - i)
                * pow(x, -i);
    }
  return result / factorial(k);
}
