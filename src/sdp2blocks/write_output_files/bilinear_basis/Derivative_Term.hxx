#pragma once

#include "../../Boost_Float.hxx"

#include <set>
#include <map>
#include <vector>

// Type to hold and compute derivatives of the form
//
// e^(-f(x)) d^k/dx^k ( e^(f(x)))

struct Derivative_Term
{
  Boost_Float constant;
  std::map<int64_t, int64_t> powers;

  Derivative_Term(const Boost_Float &Constant,
                  const std::map<int64_t, int64_t> &Powers)
      : constant(Constant), powers(Powers)
  {}
  std::set<Derivative_Term> derivative() const;
};

inline bool operator<(const Derivative_Term &a, const Derivative_Term &b)
{
  return a.powers < b.powers
           ? true
           : (a.powers > b.powers ? false : a.constant < b.constant);
}
