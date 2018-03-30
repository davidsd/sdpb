//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include <mblas.h>
#include <mlapack.h>

using Integer = mpackint;
using Real = mpf_class;

inline bool compare_abs(const Real &a, const Real &b)
{
  return abs(a) < abs(b);
}
