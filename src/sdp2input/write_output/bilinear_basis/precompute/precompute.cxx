#include "../accumulate_over_others.hxx"

#include <boost/math/special_functions/expint.hpp>

#include <vector>

Boost_Float integral(const Boost_Float &prefactor, const Boost_Float &b,
                     const Boost_Float &x, const int64_t &k);

void precompute(
  const Boost_Float &base, std::vector<Boost_Float> &sorted_poles,
  std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                        std::vector<Boost_Float>::const_iterator>>
    &equal_ranges,
  std::vector<int64_t> &lengths, std::vector<Boost_Float> &products,
  std::vector<std::vector<Boost_Float>> &integral_matrix)
{
  sorted_poles.erase(
    std::remove_if(sorted_poles.begin(), sorted_poles.end(),
                   [](const Boost_Float &a) { return a >= 0; }),
    sorted_poles.end());
  std::sort(sorted_poles.begin(), sorted_poles.end());

  Boost_Float tolerance(
    pow(Boost_Float(10.0), -(El::gmp::Precision() * std::log10(2.0)) / 2));
  for(auto pole(sorted_poles.begin()); pole != sorted_poles.end();)
    {
      const Boost_Float &p(*pole);
      equal_ranges.push_back(
        std::equal_range(pole, sorted_poles.end(), p,
                         [&](const Boost_Float &p, const Boost_Float &q) {
                           return (p < q) && !(abs(p - q) < tolerance);
                         }));
      auto &equal_range(equal_ranges.back());
      lengths.push_back(std::distance(equal_range.first, equal_range.second));
      int64_t l(lengths.back());

      products.push_back(
        1
        / accumulate_over_others(
          sorted_poles, equal_range, Boost_Float(1),
          [&](const Boost_Float &product, const Boost_Float &q) {
            return product * (p - q);
          }));

      Boost_Float integral_sum(0);
      Boost_Float integral_prefactor(-boost::math::expint(-p * log(base))
                                     * pow(base, p));

      integral_matrix.emplace_back();
      auto &integrals(integral_matrix.back());
      for(int64_t k = 0; k < l; ++k)
        {
          integrals.push_back(
            integral(integral_prefactor, base, p, l - k - 1));
        }

      std::advance(pole, l);
    }
}
