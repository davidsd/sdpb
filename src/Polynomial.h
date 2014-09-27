#ifndef SDP_BOOTSTRAP_POLYNOMIAL_H_
#define SDP_BOOTSTRAP_POLYNOMIAL_H_

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

class PolynomialVectorMatrix {
public:
  int rows;
  int cols;
  vector<vector<Polynomial> > elements;

  inline const vector<Polynomial>& elt(const int r, const int c) const {
    return elements[r + c*rows];
  }

  inline vector<Polynomial>& elt(const int r, const int c) {
    return elements[r + c*rows];
  }

  int degree() const {
    int d = 0;
    for (vector<vector<Polynomial> >::const_iterator e = elements.begin(); e != elements.end(); e++)
      for (vector<Polynomial>::const_iterator p = e->begin(); p != e->end(); p++)
        d = max(p->degree(), d);
    return d;
  }

};


#endif  // SDP_BOOTSTRAP_POLYNOMIAL_H_
