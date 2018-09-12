//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include <El.hpp>

#include <cassert>
#include <iostream>
#include <vector>

// A univariate polynomial
//
//   p(x) = a_0 + a_1 x + a_2 x^2 + ... + a_n x^n
//
class Polynomial
{
public:
  // Coefficients {a_0, a_1, ..., a_n} in increasing order of degree
  std::vector<El::BigFloat> coefficients;

  // The zero polynomial
  Polynomial() : coefficients(1, 0) {}

  // Degree of p(x)
  int degree() const { return coefficients.size() - 1; };

  // Evaluate p(x) for some x using horner's method
  El::BigFloat operator()(const El::BigFloat &x) const
  {
    assert(!coefficients.empty());
    auto coefficient(coefficients.rbegin());
    El::BigFloat result(*coefficient);
    ++coefficient;
    for(; coefficient != coefficients.rend(); ++coefficient)
      {
        result *= x;
        result += *coefficient;
      }
    return result;
  }

  // Print p(x), for debugging purposes
  friend std::ostream &operator<<(std::ostream &os, const Polynomial &p)
  {
    for(int i = p.degree(); i >= 0; i--)
      {
        os << p.coefficients[i];
        if(i > 1)
          {
            os << "x^" << i << " + ";
          }
        else if(i == 1)
          {
            os << "x + ";
          }
      }
    return os;
  }
};
