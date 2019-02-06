#include "../../../Polynomial.hxx"
#include "../../Damped_Rational.hxx"

Boost_Float
bilinear_form(const Damped_Rational &damped_rational, const int64_t &m);

std::vector<Polynomial> bilinear_basis(const Damped_Rational &damped_rational,
                                       const size_t &half_max_degree)
{
  std::vector<Polynomial> result;

  Polynomial polynomial(half_max_degree, 1);

  bilinear_form(damped_rational,2*half_max_degree);
  
  return result;
}
