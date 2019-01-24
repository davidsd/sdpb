#pragma once

#include <El.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <vector>

using Boost_Float = boost::multiprecision::mpfr_float;
struct Damped_Rational
{
  Boost_Float constant, base;
  std::vector<Boost_Float> poles;
  std::string variable;
};

inline void swap(Damped_Rational &a, Damped_Rational &b)
{
  using namespace std;
  swap(a.constant, b.constant);
  swap(a.base, b.base);
  swap(a.poles, b.poles);
  swap(a.variable, b.variable);
}
