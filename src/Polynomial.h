//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_POLYNOMIAL_H_
#define SDPB_POLYNOMIAL_H_

#include "types.h"
#include "util.h"
#include "Vector.h"

class Polynomial {
public:
  Vector coefficients;

  Polynomial(): coefficients(1, 0) {}

  int degree() const {
    return coefficients.size() - 1;
  };

  Real operator()(const Real &x) const {
    int deg = degree();
    Real y = coefficients[deg];
    for (int i = deg - 1; i >= 0; i--) {
      y *= x;
      y += coefficients[i];
    }
    return y;
  }

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
