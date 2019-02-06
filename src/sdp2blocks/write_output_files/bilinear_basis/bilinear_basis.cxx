#include "../../../Polynomial.hxx"
#include "../../Damped_Rational.hxx"

#include "Derivative_Term.hxx"

#include <set>

std::set<Derivative_Term> dExp(const int64_t &k);

std::vector<Polynomial> bilinear_basis(const Damped_Rational &damped_rational,
                                       const size_t &half_max_degree)
{
  std::vector<Polynomial> result;

  Polynomial polynomial(half_max_degree, 1);

  std::cout << "About to Exp\n";
  for(auto &d : dExp(10))
    {
      std::cout << d << "\n";
    }
  std::cout << std::flush;

  return result;
}
