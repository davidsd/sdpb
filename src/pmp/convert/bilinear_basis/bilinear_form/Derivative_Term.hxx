#pragma once

#include "sdpb_util/Boost_Float.hxx"

#include <map>
#include <vector>
#include <iostream>

// Type to hold and compute derivatives of the form
//
// e^(-f(x)) d^k/dx^k ( e^(f(x)))

// A single term is of the form
//
// (constant) (d(^a0)f)^b0 (d(^a1)f)^b1 (d(^a2)f)^b2 ... 
//
// The complete derivative is a collection of these terms

struct Derivative_Term
{
  Boost_Float constant;
  // powers: key is the n'th derivative, value is the exponent.  So an
  // entry [5,3] is the 5th derivative cubed.
  std::map<int64_t, int64_t> powers;

  Derivative_Term(const Boost_Float &Constant,
                  const std::map<int64_t, int64_t> &Powers)
      : constant(Constant), powers(Powers)
  {}
};

inline bool operator<(const Derivative_Term &a, const Derivative_Term &b)
{
  return std::tie(a.powers, a.constant) < std::tie(b.powers, b.constant);
}

inline std::ostream & operator<<(std::ostream &os, const Derivative_Term &term)
{
  os << term.constant;
  for(auto &power: term.powers)
    {
      os << " * d(" << power.first << ",f)";
      if(power.second!=1)
        {
          os << "^" << power.second;
        }
    }
  return os;
}
