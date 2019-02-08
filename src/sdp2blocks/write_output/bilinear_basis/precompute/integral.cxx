#include "../factorial.hxx"

// integral(0) = base^x tgamma(0,x*log(base))
//             = -(base^x) expint(0, -x log(base))
// integral(k) = d(integral(k-1),x)/k
//
// In closed form
//
// integral(k) = -(base^x expint(0, -x log(base)) * (log(base))^k
//                + Sum(i={1..=k}) {(-1)^k * (i-1)! (-log(base))^(k-i) / x^i})/k!

Boost_Float integral(const Boost_Float &prefactor, const Boost_Float &base,
                     const Boost_Float &x, const int64_t &k)
{
  Boost_Float result(prefactor * pow(log(base), k));
  for(int64_t i = 1; i <= k; ++i)
    {
      result += (k % 2 == 0 ? 1 : -1) * factorial(i - 1)
                * pow(-log(base), k - i) * pow(x, -i);
    }
  return result / factorial(k);
}
