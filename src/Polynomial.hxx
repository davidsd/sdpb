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

// Convenience functions to avoid copies
inline void swap(std::vector<Polynomial> &polynomials,
                 std::vector<std::vector<El::BigFloat>> &elements_vector)
{
  for(auto &elements : elements_vector)
    {
      polynomials.emplace_back();
      std::swap(polynomials.back().coefficients, elements);
    }
}

inline void swap(
  std::vector<std::vector<Polynomial>> &polynomials_vector,
  std::vector<std::vector<std::vector<El::BigFloat>>> &elements_vector_vector)
{
  for(auto &elements_vector : elements_vector_vector)
    {
      polynomials_vector.emplace_back();
      swap(polynomials_vector.back(), elements_vector);
    }
}
