//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_VECTOR_H_
#define SDPB_VECTOR_H_

#include <algorithm>
#include <vector>
#include <assert.h>
#include "types.h"
#include "util.h"

using std::vector;

typedef vector<Real> Vector;

inline Real maxAbsVector(const Vector &v) {
  return abs(*std::max_element(v.begin(), v.end(), compareAbs));
}  

inline void fillVector(Vector &v, const Real &a) {
  std::fill(v.begin(), v.end(), a);
}

inline void scaleVector(Vector &v, const Real &a) {
  for (unsigned int i = 0; i < v.size(); i++)
    v[i] *= a;
}

inline void addVector(Vector &v, const Vector &u) {
  assert(v.size() == u.size());

  for (unsigned int i = 0; i < v.size(); i++)
    v[i] += u[i];
}

inline Real dotProduct(const Vector &u, const Vector v) {
  Real result = 0;
  for (unsigned int i = 0; i < u.size(); i++)
    result += u[i]*v[i];
  return result;
}

inline Vector multiplyVectors(const Vector &u, const Vector &v) {
  Vector w(u.size());
  for (unsigned int i = 0; i < w.size(); i++)
    w[i] = u[i]*v[i];
  return w;
}

#endif  // SDPB_VECTOR_H_
