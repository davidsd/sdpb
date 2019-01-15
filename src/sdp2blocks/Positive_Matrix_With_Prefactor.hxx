#pragma once

#include "../Polynomial.hxx"

struct Positive_Matrix_With_Prefactor
{
  Damped_Rational damped_rational;
  std::vector<Polynomial> polynomials;
};

inline void
swap(Positive_Matrix_With_Prefactor &a, Positive_Matrix_With_Prefactor &b)
{
  swap(a.damped_rational, b.damped_rational);
  swap(a.polynomials, b.polynomials);
}
