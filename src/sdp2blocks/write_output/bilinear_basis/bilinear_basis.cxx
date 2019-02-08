#include "../../Damped_Rational.hxx"
#include "../../../Polynomial.hxx"
#include "../../../Timers.hxx"

void precompute(
  const Boost_Float &base, std::vector<Boost_Float> &sorted_poles,
  std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                        std::vector<Boost_Float>::const_iterator>>
    &equal_ranges,
  std::vector<int64_t> &lengths, std::vector<Boost_Float> &products,
  std::vector<std::vector<Boost_Float>> &integral_matrix);

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
  std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                        std::vector<Boost_Float>::const_iterator>>
    equal_ranges;
  std::vector<int64_t> lengths;
  std::vector<Boost_Float> products;
  std::vector<std::vector<Boost_Float>> integral_matrix;
  precompute(damped_rational.base, sorted_poles, equal_ranges, lengths,
             products, integral_matrix);

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
