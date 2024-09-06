#pragma once

#include "Polynomial.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"

#include <El.hpp>
#include <vector>

struct Polynomial_Power
{
  Boost_Polynomial polynomial{};
  Boost_Float power = 1;

  [[nodiscard]] Boost_Float evaluate(const Boost_Float &x) const
  {
    // TODO too many conversions! Make everything Boost_Float?
    return pow(polynomial(x), power);
  }

  void clear()
  {
    polynomial ={};
    power = 1;
  }
};

struct Polynomial_Power_Product
{
  std::vector<Polynomial_Power> terms{};

  [[nodiscard]] Boost_Float evaluate(const Boost_Float &x) const
  {
    Boost_Float result = 1;
    for(const auto &term : terms)
      {
        result *= term.evaluate(x);
      }
    return result;
  }

  [[nodiscard]] El::BigFloat evaluate(const El::BigFloat &x) const
  {
    if(terms.empty())
      return 1;
    return to_BigFloat(evaluate(to_Boost_Float(x)));
  }
};
