#include "accumulate_over_others.hxx"
#include "../../Damped_Rational.hxx"
#include "../../../Polynomial.hxx"
#include "../../../Timers.hxx"

#include <boost/math/special_functions/expint.hpp>

Boost_Float integral(const Boost_Float &prefactor, const Boost_Float &b,
                     const Boost_Float &x, const int64_t &k);

Boost_Float bilinear_form(
  const Damped_Rational &damped_rational,
  const std::vector<Boost_Float> &sorted_poles,
  const std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                              std::vector<Boost_Float>::const_iterator>>
    &equal_ranges,
  const std::vector<int64_t> &lengths,
  const std::vector<Boost_Float> &products,
  const std::vector<std::vector<Boost_Float>> &integral_matrix,
  const int64_t &m);

std::vector<Polynomial>
bilinear_basis(const Damped_Rational &damped_rational,
               const size_t &half_max_degree, const std::string &timer_prefix,
               Timers &timers)
{
  std::vector<Polynomial> result;

  Polynomial polynomial(half_max_degree, 1);

  std::vector<Boost_Float> sorted_poles(damped_rational.poles);
  sorted_poles.erase(
    std::remove_if(sorted_poles.begin(), sorted_poles.end(),
                   [](const Boost_Float &a) { return a >= 0; }),
    sorted_poles.end());
  std::sort(sorted_poles.begin(), sorted_poles.end());

  std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                        std::vector<Boost_Float>::const_iterator>>
    equal_ranges;
  std::vector<int64_t> lengths;
  std::vector<Boost_Float> products;
  std::vector<std::vector<Boost_Float>> integral_matrix;
  for(auto pole(sorted_poles.begin()); pole != sorted_poles.end();)
    {
      const Boost_Float &p(*pole);
      equal_ranges.push_back(
        std::equal_range(pole, sorted_poles.end(), p,
                         [&](const Boost_Float &p, const Boost_Float &q) {
                           return (p < q) && !(abs(p - q) < 1e-2);
                         }));
      auto &equal_range(equal_ranges.back());
      lengths.push_back(std::distance(equal_range.first, equal_range.second));
      int64_t l(lengths.back());
      
      products.push_back(1/accumulate_over_others(
        sorted_poles, equal_range, Boost_Float(1),
        [&](const Boost_Float &product, const Boost_Float &q) {
          return product * (p - q);
        }));

      Boost_Float integral_sum(0);
      Boost_Float integral_prefactor(
        -boost::math::expint(-p * log(damped_rational.base))
        * pow(damped_rational.base, p));

      integral_matrix.emplace_back();
      auto &integrals(integral_matrix.back());
      for(int64_t k = 0; k < l; ++k)
        {
          integrals.push_back(
            integral(integral_prefactor, damped_rational.base, p, l - k - 1));
        }

      std::advance(pole, l);
    }

  std::vector<Boost_Float> bilinear_table;
  for(int64_t m = 0; m <= int64_t(2 * half_max_degree); ++m)
    {
      const std::string bilinear_form_timer_name(
        timer_prefix + ".bilinear_form_" + std::to_string(m));
      auto &bilinear_form_timer(
        timers.add_and_start(bilinear_form_timer_name));

      bilinear_table.push_back(bilinear_form(damped_rational, sorted_poles,
                                             equal_ranges, lengths, products,
                                             integral_matrix, m));

      bilinear_form_timer.stop();

      if(m == 38)
        {
          std::cout.precision(Boost_Float::default_precision());
          std::cout << "m: " << m << " " << bilinear_table.back() << "\n";
        }
    }
  return result;
}
