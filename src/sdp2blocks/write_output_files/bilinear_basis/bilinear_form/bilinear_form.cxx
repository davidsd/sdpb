#include "accumulate_over_others.hxx"
#include "../../../Damped_Rational.hxx"
#include "../../../../Polynomial.hxx"

Boost_Float
integral(const Boost_Float &b, const Boost_Float &x, const int64_t &k);

Boost_Float
rest(const int64_t &m, const Boost_Float &p,
     const std::vector<Boost_Float> &sorted_poles,
     const std::pair<std::vector<Boost_Float>::const_iterator,
                     std::vector<Boost_Float>::const_iterator> &equal_range,
     const int64_t &k);

Boost_Float
bilinear_form(const Damped_Rational &damped_rational, const int64_t &m)
{
  Boost_Float result(damped_rational.constant);
  std::vector<Boost_Float> sorted_poles(damped_rational.poles);

  sorted_poles.erase(
    std::remove_if(sorted_poles.begin(), sorted_poles.end(),
                   [](const Boost_Float &a) { return a >= 0; }),
    sorted_poles.end());
  std::sort(sorted_poles.begin(), sorted_poles.end());

  Boost_Float pole_sum(0);
  for(auto pole(sorted_poles.begin()); pole != sorted_poles.end();)
    {
      Boost_Float &p(*pole);
      const std::pair<std::vector<Boost_Float>::const_iterator,
                      std::vector<Boost_Float>::const_iterator>
        equal_range(std::equal_range(pole, sorted_poles.end(), p));
      auto l(std::distance(equal_range.first, equal_range.second));

      Boost_Float product(accumulate_over_others(
        sorted_poles, equal_range, Boost_Float(1),
        [&](const Boost_Float &product, const Boost_Float &q) {
          return product * (p - q);
        }));

      Boost_Float integral_sum(0);
      for(int64_t k = 0; k < l; ++k)
        {
          integral_sum += integral(damped_rational.base, p, l - k - 1)
                          * rest(m, p, sorted_poles, equal_range, k);
        }
      pole_sum += (pow(p, m) / product) * integral_sum;
      do
        {
          ++pole;
        }
      while(pole != sorted_poles.end() && abs(p - *pole) < 1.0e-2);
    }
  std::cout << "b: " << damped_rational.base << "\n"
            << "m: " << m << "\n"
            << "pole_sum: " << pole_sum << "\n"
            << std::flush;
}
