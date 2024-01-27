//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <iostream>
#include <vector>
#include <boost/math/tools/polynomial.hpp>

// FIXME: Use boost::math::tools::polynomial instead

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
  Polynomial(const size_t &size, const El::BigFloat &default_element)
      : coefficients(size, default_element)
  {}

  // Degree of p(x)
  int64_t degree() const { return coefficients.size() - 1; };

  // Evaluate p(x) for some x using horner's method
  El::BigFloat operator()(const El::BigFloat &x) const
  {
    ASSERT(!coefficients.empty());
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

  friend bool operator==(const Polynomial &lhs, const Polynomial &rhs)
  {
    return lhs.coefficients == rhs.coefficients;
  }
  friend bool operator!=(const Polynomial &lhs, const Polynomial &rhs)
  {
    return !(lhs == rhs);
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

inline Polynomial operator/(const Polynomial &a, const El::BigFloat &b)
{
  Polynomial result(0, 0);
  result.coefficients.reserve(a.degree() + 1);
  const El::BigFloat inverse(1 / b);
  for(auto &coefficient : a.coefficients)
    {
      result.coefficients.emplace_back(coefficient * inverse);
    }
  return result;
}

struct Polynomial_Vector : std::vector<Polynomial>
{
  template <class... Args>
  explicit Polynomial_Vector(Args &&...args)
      : std::vector<Polynomial>(std::forward<Args>(args)...)
  {}

  // Need to implement this method to use El::Matrix<Polynomial_Vector>:
  // for some configurations, El::MemZero() is called when constructing a matrix.
  // It calls element.Zero() for each matrix element.
  void Zero() { clear(); }
};

// Convenience functions to avoid copies
inline void swap(Polynomial_Vector &polynomials,
                 std::vector<std::vector<El::BigFloat>> &elements_vector)
{
  for(auto &elements : elements_vector)
    {
      polynomials.emplace_back();
      std::swap(polynomials.back().coefficients, elements);
    }
}

inline void swap(
  std::vector<Polynomial_Vector> &polynomials_vector,
  std::vector<std::vector<std::vector<El::BigFloat>>> &elements_vector_vector)
{
  for(auto &elements_vector : elements_vector_vector)
    {
      polynomials_vector.emplace_back();
      swap(polynomials_vector.back(), elements_vector);
    }
}

inline void
swap(std::vector<std::vector<Polynomial_Vector>> &polynomials_vector_vector,
     std::vector<std::vector<std::vector<std::vector<El::BigFloat>>>>
       &elements_vector_vector_vector)
{
  for(auto &elements_vector_vector : elements_vector_vector_vector)
    {
      polynomials_vector_vector.emplace_back();
      swap(polynomials_vector_vector.back(), elements_vector_vector);
    }
}
