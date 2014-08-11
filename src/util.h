#ifndef SDP_BOOTSTRAP_UTIL_H_
#define SDP_BOOTSTRAP_UTIL_H_

#include <algorithm>
#include <iostream>
#include <ostream>
#include <vector>
#include "types.h"

using std::ostream;
using std::vector;

inline bool compareAbs(const Real &a, const Real &b) {
  return abs(a) < abs(b);
}

template <class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  os << "{";
  int last = v.size() - 1;
  for (int i = 0; i < last; i++)
    os << v[i] << ", ";
  if (last >= 0)
    os << v[last];
  os << "}";
  return os;
}

template<class Iter, class T>
Iter binaryFind(Iter begin, Iter end, T val)
{
  // Finds the lower bound in at most log(last - first) + 1 comparisons
  Iter i = std::lower_bound(begin, end, val);

  if (i != end && !(val < *i))
    return i; // found
  else
    return end; // not found
}

#endif  // SDP_BOOTSTRAP_UTIL_H_
