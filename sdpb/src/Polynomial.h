//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_POLYNOMIAL_H_
#define SDPB_POLYNOMIAL_H_

#include "types.h"
#include "Vector.h"

// A univariate polynomial
// 
//   p(x) = a_0 + a_1 x + a_2 x^2 + ... + a_n x^n
//
class Polynomial {
 public:
  // Coefficients {a_0, a_1, ..., a_n} in increasing order of degree
  Vector coefficients;

  // The zero polynomial
  Polynomial(): coefficients(1, 0) {}

  // Degree of p(x)
  int degree() const {
    return coefficients.size() - 1;
  };

  // Evaluate p(x) for some x using horner's method
  Real operator()(const Real &x) const {
    int deg = degree();
    Real result = coefficients[deg];
    for (int i = deg - 1; i >= 0; i--) {
      result *= x;
      result += coefficients[i];
    }
    return result;
  }

  // Print p(x), for debugging purposes
  friend ostream& operator<<(ostream& os, const Polynomial& p) {
    for (int i = p.degree(); i >= 0; i--) {
      os << p.coefficients[i];
      if (i > 1)
        os << "x^" << i << " + ";
      else if (i == 1)
        os << "x + ";
    }
    return os;
  }
};

#endif  // SDPB_POLYNOMIAL_H_
