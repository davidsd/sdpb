//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_UTIL_H_
#define SDPB_UTIL_H_

#include <algorithm>
#include <iostream>
#include <ostream>
#include <vector>
#include "types.h"

using std::ostream;
using std::vector;

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

#endif  // SDPB_UTIL_H_
