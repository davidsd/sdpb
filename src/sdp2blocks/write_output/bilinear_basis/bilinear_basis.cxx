#include "../../../Polynomial.hxx"
#include "../../Damped_Rational.hxx"

Boost_Float
bilinear_form(const Damped_Rational &damped_rational, const int64_t &m);

std::vector<Polynomial> bilinear_basis(const Damped_Rational &damped_rational,
                                       const size_t &half_max_degree)
{
  std::vector<Polynomial> result;

  Polynomial polynomial(half_max_degree, 1);

  std::vector<Boost_Float> bilinear_table;
  for(int64_t m = 0; m <= 2 * half_max_degree; ++m)
    {
      bilinear_table.push_back(bilinear_form(damped_rational, m));
    }
  
  return result;
}
