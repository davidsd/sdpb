#pragma once

#include "Damped_Rational.hxx"
#include "sdpb_util/Polynomial.hxx"

struct Positive_Matrix_With_Prefactor
{
  Damped_Rational damped_rational;
  // TODO: This is actually a symmetric matrix (dim x dim) of a vector
  // of polynomials.  So we should use a proper symmetric matrix class.
  std::vector<std::vector<std::vector<Polynomial>>> polynomials;

  // Prohibit accidential copying to avoid excessive memory usage,
  // allow only moving

  Positive_Matrix_With_Prefactor() = default;

  Positive_Matrix_With_Prefactor(const Positive_Matrix_With_Prefactor &other)
    = delete;
  Positive_Matrix_With_Prefactor(
    Positive_Matrix_With_Prefactor &&other) noexcept
    = default;

  Positive_Matrix_With_Prefactor &
  operator=(const Positive_Matrix_With_Prefactor &other)
    = delete;

  Positive_Matrix_With_Prefactor &
  operator=(Positive_Matrix_With_Prefactor &&other) noexcept
    = default;
};

inline void
swap(Positive_Matrix_With_Prefactor &a, Positive_Matrix_With_Prefactor &b)
{
  swap(a.damped_rational, b.damped_rational);
  swap(a.polynomials, b.polynomials);
}
